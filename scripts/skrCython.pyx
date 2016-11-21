import errno

cpdef getDistance(gene_name, start, end, gene_strand, dic, mode):
    pos_ranges = dic[gene_name]
    if gene_strand == "-": pos_ranges = pos_ranges[::-1]
    num_regions = len(pos_ranges)
    if mode == 5:
        pos = start + 15 if gene_strand == "+" else end - 15
        return offsetHelper(num_regions, pos_ranges, gene_strand, pos)
    elif mode == 3:
        pos = end - 12 if gene_strand == "+" else start + 12
        return offsetHelper(num_regions, pos_ranges, gene_strand, pos)
    else:
        exit(errno.EINVAL)

cpdef offsetHelper(numRegions, posRanges, geneStrand, pos):
        i = 0
        while i < numRegions:
            left, right, phase = posRanges[i]
            if left <= pos <= right:
                up = left if geneStrand == "+" else right
                offset = (abs(pos - up) + phase) % 3
                return offset
            else:
                i += 1
        return -1