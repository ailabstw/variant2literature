"""detect tables in pdf
"""
import logging
import os
import threading
import io
import time
import random

import cv2
import GPUtil
import rpyc
from rpyc.utils.helpers import classpartial
import numpy as np
import torch
from torch import multiprocessing
from torch.autograd import Variable

import _init_paths
from model.nms.nms_wrapper import nms
from model.faster_rcnn.vgg16 import vgg16
from model.faster_rcnn.resnet import resnet
from model.utils.config import cfg
from model.utils.blob import im_list_to_blob
from model.rpn.bbox_transform import clip_boxes
from model.rpn.bbox_transform import bbox_transform_inv

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def load_np(data):
    """load dumped numpy array
    """
    return np.load(io.BytesIO(data))['data']


def dump_np(data):
    """dump numpy array
    """
    f = io.BytesIO()
    np.savez(f, data=data)
    f.seek(0)
    return f.read()


def get_image_blob(im):
    """Converts an image into a network input.
    Arguments:
        im (ndarray): a color image in BGR order
    Returns:
        blob (ndarray): a data blob holding an image pyramid
        im_scale_factors (list): list of image scales (relative to im) used
                                 in the image pyramid
    """
    im_orig = im.astype(np.float32, copy=True)
    im_orig -= cfg.PIXEL_MEANS

    im_shape = im_orig.shape
    im_size_min = np.min(im_shape[0:2])
    im_size_max = np.max(im_shape[0:2])

    processed_ims = []
    im_scale_factors = []

    for target_size in cfg.TEST.SCALES:
        im_scale = float(target_size) / float(im_size_min)
    # Prevent the biggest axis from being more than MAX_SIZE
    if np.round(im_scale * im_size_max) > cfg.TEST.MAX_SIZE:
        im_scale = float(cfg.TEST.MAX_SIZE) / float(im_size_max)
    im = cv2.resize(im_orig, None, None, fx=im_scale, fy=im_scale,
                    interpolation=cv2.INTER_LINEAR)
    im_scale_factors.append(im_scale)
    processed_ims.append(im)

    # Create a blob to hold the input images
    blob = im_list_to_blob(processed_ims)
    return blob, np.array(im_scale_factors)


class FasterRCNN:
    """Faster RCNN pdf table detector
    """

    def load_model(self):
        checkpoint_path = '/app/models/faster_rcnn.pth'
        self.classes = np.array(['__background__', 'table'])
        self.model = vgg16(self.classes, pretrained=False, class_agnostic=False)
        self.model.create_architecture()

        checkpoint = torch.load(checkpoint_path)
        self.model.load_state_dict(checkpoint['model'])

        if 'pooling_mode' in checkpoint.keys():
            cfg.POOLING_MODE = checkpoint['pooling_mode']

        self.im_data = torch.FloatTensor(1)
        self.im_info = torch.FloatTensor(1)
        self.num_boxes = torch.LongTensor(1)
        self.gt_boxes = torch.FloatTensor(1)

        self.im_data = self.im_data.cuda()
        self.im_info = self.im_info.cuda()
        self.num_boxes = self.num_boxes.cuda()
        self.gt_boxes = self.gt_boxes.cuda()
        self.model.cuda()

        with torch.no_grad():
            self.im_data = Variable(self.im_data)
            self.im_info = Variable(self.im_info)
            self.num_boxes = Variable(self.num_boxes)
            self.gt_boxes = Variable(self.gt_boxes)
        self.model.eval()

    # def detect(self, blobs, im_scales):
    def detect(self, img_data):
        np_arr = np.fromstring(img_data, np.uint8)
        image = cv2.imdecode(np_arr, cv2.IMREAD_COLOR)
        blobs, im_scales = get_image_blob(image)

        im_blob = blobs
        im_info_np = np.array([[im_blob.shape[1], im_blob.shape[2], im_scales[0]]], dtype=np.float32)
        im_data_pt = torch.from_numpy(im_blob)
        im_data_pt = im_data_pt.permute(0, 3, 1, 2)
        im_info_pt = torch.from_numpy(im_info_np)

        self.im_data.data.resize_(im_data_pt.size()).copy_(im_data_pt)
        self.im_info.data.resize_(im_info_pt.size()).copy_(im_info_pt)
        self.gt_boxes.data.resize_(1, 1, 5).zero_()
        self.num_boxes.data.resize_(1).zero_()

        rois, cls_prob, bbox_pred, *_ = \
            self.model(self.im_data, self.im_info, self.gt_boxes, self.num_boxes)

        scores = cls_prob.data
        boxes = rois.data[:, :, 1:5]

        box_deltas = bbox_pred.data
        box_deltas = box_deltas.view(-1, 4) * torch.FloatTensor(cfg.TRAIN.BBOX_NORMALIZE_STDS).cuda() \
                   + torch.FloatTensor(cfg.TRAIN.BBOX_NORMALIZE_MEANS).cuda()
        box_deltas = box_deltas.view(1, -1, 4 * len(self.classes))

        pred_boxes = bbox_transform_inv(boxes, box_deltas, 1)
        pred_boxes = clip_boxes(pred_boxes, self.im_info.data, 1)

        pred_boxes /= im_scales[0]

        scores = scores.squeeze()
        pred_boxes = pred_boxes.squeeze()

        thresh = 0.7
        inds = torch.nonzero(scores[:, 1] > thresh).view(-1)
        cls_dets = []
        if inds.numel() > 0:
            cls_scores = scores[:, 1][inds]
            _, order = torch.sort(cls_scores, 0, True)
            cls_boxes = pred_boxes[inds][:, 4:8]

            cls_dets = torch.cat((cls_boxes, cls_scores.unsqueeze(1)), 1)
            cls_dets = cls_dets[order]
            keep = nms(cls_dets, cfg.TEST.NMS, force_cpu=False)
            cls_dets = cls_dets[keep.view(-1).long()]
            cls_dets = cls_dets.cpu().numpy()
        cls_dets = [cls_det[:4] for cls_det in cls_dets]
        return cls_dets


def worker(que):
    faster_rcnn = FasterRCNN()
    faster_rcnn.load_model()
    logger.info('init OK')

    while True:
        msg = que.get()
        img_data, pipe = msg
        results = faster_rcnn.detect(img_data)
        pipe.send(dump_np(results))


class TableDetector(rpyc.Service):
    """rpyc service
    """

    def __init__(self, que):
        super(TableDetector, self).__init__()
        self.que = que
        self.pipe = multiprocessing.Pipe()

    def exposed_detect(self, img_data):
        """detect tables in image
        """
        t = time.time()
        self.que.put((img_data, self.pipe[1]))
        ret = self.pipe[0].recv()
        dt = time.time() - t
        logger.info(f'tables detected. {dt:.3f} secs')
        return ret


def main():
    """main
    """
    multiprocessing.set_start_method('spawn', force=True)
    que = multiprocessing.Queue()

    n_process = int(os.environ.get('NUM_TABLE_DETECTORS', '1'))
    for i in range(n_process):
        p = multiprocessing.Process(target=worker, args=(que,), daemon=True)
        p.start()

    service = classpartial(TableDetector, que=que)
    t = rpyc.utils.server.ThreadedServer(service, port=18861)
    t.start()


if __name__ == '__main__':
    main()
