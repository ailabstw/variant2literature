import re


class Node:
    def __init__(self, is_end=False, gene_id=''):
        self.nxt = dict()
        self.is_end = is_end
        self.gene_id = gene_id


class GenePrefixTree:
    def __init__(self, path='/app/models/PT_CTDGene.txt'):
        self.root = Node()
        path_dict = {(): self.root}

        with open(path) as f:
            f.read(1)
            for line in f:
                fields = line.strip('\n').split('\t')
                path, token = fields[:2]
                is_end = len(fields) >= 3
                gene_id = fields[2] if is_end else ''

                idxes = tuple(path.split('-'))
                parent = path_dict.get(idxes[:-1])
                path_dict[idxes] = Node(is_end, gene_id)
                parent.nxt[token] = path_dict[idxes]

    def remove_punkt(self, tokens, offsets):
        new_tokens, new_offsets = [], []
        for token, offsets in zip(tokens, offsets):
            if not re.search(r'\W', token):
                new_tokens.append(token)
                new_offsets.append(offsets)
        return new_tokens, new_offsets

    def search_tokens(self, tokens, offsets):
        tokens, offsets = self.remove_punkt(tokens, offsets)
        i, genes = 0, []
        while i < len(tokens):
            j, node = 0, self.root
            start, end = i, None
            while node and i + j < len(tokens):
                node = node.nxt.get(tokens[i + j].lower(), None)
                if node and node.is_end:
                    end = i + j
                j += 1

            if end:
                genes.append((offsets[start][0], offsets[end][1]))
                i = end
            i += 1
        return genes

    def search(self, tokens):
        node = self.root
        for token in tokens:
            node = node.nxt.get(token, None)
            if not node:
                return False
        return node.is_end


if __name__ == '__main__':
    gene_pt = GenePrefixTree()
    print(gene_pt.search(['myo', '15', 'a']))
