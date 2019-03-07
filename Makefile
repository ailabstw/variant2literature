CURRENT_DIR = $(shell pwd)
DOCKER_HOST = $(shell ip -4 addr show docker0 | grep -Po 'inet \K[\d.]+')

MYSQL_PORT=55688
MYSQL_ROOT_PASSWORD=s8fjYJd92oP
MYSQL_VOLUME=${CURRENT_DIR}/mysql_data

CUDA_VISIBLE_DEVICES=0
NUM_PROCESSES=1
NUM_TABLE_DETECTORS=1

IMAGE_NAME=variant2literature
CONTAINER_NAME=v2l
MYSQL_NAME=v2l_mysql

build:
	docker build -t ${IMAGE_NAME} .

compile:
	nvidia-docker run --rm --name ${CONTAINER_NAME} \
		-v ${CURRENT_DIR}:/app \
		${IMAGE_NAME} \
		bash -c "cd table_detector/lib && bash make.sh"

run:
	nvidia-docker run -d --name ${CONTAINER_NAME} \
		-v ${CURRENT_DIR}:/app \
		-e MYSQL_HOST=${DOCKER_HOST} \
		-e MYSQL_PORT=${MYSQL_PORT} \
		-e MYSQL_ROOT_PASSWORD=${MYSQL_ROOT_PASSWORD} \
		-e CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES} \
		-e NUM_TABLE_DETECTORS=${NUM_TABLE_DETECTORS} \
		-e LOAD_BALANCER_HOST='localhost' \
		${IMAGE_NAME} \
		bash -c "cd table_detector && python table_detector.py"

run-db:
	docker run -d --name ${MYSQL_NAME} \
		-v ${MYSQL_VOLUME}:/var/lib/mysql \
		-p ${MYSQL_PORT}:3306 \
		-e MYSQL_ROOT_PASSWORD=${MYSQL_ROOT_PASSWORD} \
		mariadb:10.3.9-bionic

load-db:
	docker exec -it ${CONTAINER_NAME} \
		bash -c "cd mysqldb && python models.py"

bash:
	docker exec -it ${CONTAINER_NAME} bash

index:
	docker exec -it ${CONTAINER_NAME} python main.py --n-process ${NUM_PROCESSES}

query:
	docker exec -it ${CONTAINER_NAME} python query.py ${OUTPUT_FILE}

rm:
	docker stop ${CONTAINER_NAME}
	docker rm ${CONTAINER_NAME}

rm-db:
	docker stop ${MYSQL_NAME}
	docker rm ${MYSQL_NAME}

truncate:
	docker exec -it ${CONTAINER_NAME} \
		bash -c "cd mysqldb && python truncate.py"
