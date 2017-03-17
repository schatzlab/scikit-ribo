# import dask.dataframe as dd
sys.path.insert(1, [i for i in sys.path if 'statsmodels' in i and '0.8.0' in i][0])

def filterRegion(self):
    # create bedtool, filter and sort
    self.startCodons = pbt.BedTool(self.startCodons).filter(lambda x: x.chrom != 'chrM')
    self.startCodons = self.startCodons.sort()
    self.orf = pbt.BedTool(self.orf).filter(lambda x: x.chrom != 'chrM')
    self.orf = self.orf.sort()
    # find overlapping regions
    distinctStartCodons = self.startCodons.merge(c="1,2", o="count").filter(lambda x: int(x[4]) == 1)
    distinctOrfs = self.orf.merge(c="1,2", o="count").filter(lambda x: int(x[4]) == 1)
    distinctStartCodons = distinctStartCodons.intersect(distinctOrfs, wa=True, sorted=True)
    # filter start codon
    startCodonHash = set([(i[0], i[1], i[2]) for i in distinctStartCodons])
    self.startCodons = self.startCodons.filter(lambda x: (x[0], x[1], x[2]) in startCodonHash)
    # filter orf
    # distinctOrfs = self.orf.merge(c="1,2", o="count").filter(lambda x: int(x[4]) == 1)
    # orfHash = set([(i[0], i[1], i[2]) for i in distinctOrfs])
    # self.orf = self.orf.filter(lambda x: (x[0], x[1], x[2]) in orfHash)

## TODO: add a function to see whether there are too many soft-clipped alignment
def filterBam(self):
    # create a template/header from the input bam file
    inBamHdl = pysam.AlignmentFile(self.inBam, "rb")
    outBamHdl = pysam.AlignmentFile(self.outBam, "wb", template=inBamHdl)
    ## read a bam file and extract info
    cigar_to_exclude = set([1,2,3,4,5]) #set(['I','D','S','H'])
    for read in inBamHdl.fetch():
        cigars = set([c[0] for c in read.cigartuples])
        # start_nt, end_nt = read.query_sequence[:2], read.query_sequence[-2:][::-1]
        # edges = set(list(start_nt + end_nt))
        # filter the bam file
        if read.mapping_quality > self.mapq and \
        self.minRL <= read.query_length <= self.maxRL and \
        len(cigars.intersection(cigar_to_exclude)) == 0 and \
        read.reference_id != 'chrM': # and 'N' not in edges:
            read.mapping_quality = read.query_length
            outBamHdl.write(read)
            # self.reads.append([read.query_name, read.query_length]) # start_nt, end_nt])
    inBamHdl.close()
    outBamHdl.close()
    print("[status]\tFinished filtering the bam file", flush=True)
    # save the bedtool to local
    self.bedtool = pbt.BedTool(self.outBam)
    self.bedtool = self.bedtool.bam_to_bed(bed12=True)
    self.bedtool.saveas(self.outBam + '.bed')


def posIndex(self):
    ## create a dict to store the position read-frame and index info
    self.posOffsets = []
    self.negOffsets = []
    with open(self.posRanges, 'r') as f:
        next(f)
        for line in f:
            gene, chr, strand, ranges = line.rstrip("\n").split("\t")
            boxes = [(int(i[0]), int(i[1]), int(i[2])) for i in [j.split(",") for j in ranges.split("|")]]
            if strand == "+":
                self.posOffsets.extend([[chr, pos, (abs(pos - (box[0] - 15)) + box[2]) % 3] for box in boxes for pos in
                                        range(box[0] - 15, box[1] + 12)])
            else:
                boxes = boxes[::-1]  # flip the order
                self.negOffsets.extend([[chr, pos, (abs(pos - (box[1] + 15)) + box[2]) % 3] for box in boxes for pos in
                                        range(box[1] + 15, box[0] - 12, -1)])
    ## convert to dataframe
    self.posOffsets = pd.DataFrame(self.posOffsets, columns=["chrom", "pos", "offset"]).drop_duplicates(
        subset=["chrom", "pos"])
    self.negOffsets = pd.DataFrame(self.negOffsets, columns=["chrom", "pos", "offset"]).drop_duplicates(
        subset=["chrom", "pos"])
    # self.offsets.to_csv(path_or_buf='offsets.testing.txt', sep='\t', header=True, index=False)
    ## save for dask
    # self.offsets = dd.from_pandas(self.offsets, npartitions=2)

def makeTrainingData(self):
    # intersect with start codons
    self.bedtool = pbt.BedTool(self.outBam + '.bed')
    trainingData = self.bedtool.intersect(self.startCodons, wa=True, wb=True, sorted=True)
    time = str(datetime.now())
    print("[status]\tFinished intersecting the bedtool with start codons:", time, flush=True)
    # convert bedtool to df
    trainingDf = trainingData.to_dataframe(names=['chrom', 'start', 'end', 'name', 'read_length', 'strand',
                                                  'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes',
                                                  'blockStarts', 'sc_chrom', 'sc_start', 'sc_end', 'gene_name',
                                                  'sc_score', 'gene_strand'])
    ### a-site
    if not self.RelE:
        trainingDf['asite'] = np.where(trainingDf['gene_strand'] == '+',
                                       trainingDf['sc_start'] - trainingDf['start'] + 3,
                                       trainingDf['end'] - trainingDf['sc_end'] + 3)
    else:
        trainingDf['asite'] = np.where(trainingDf['gene_strand'] == '+',
                                       trainingDf['end'] - trainingDf['sc_start'] - 3,
                                       trainingDf['sc_end'] - trainingDf['start'] - 3)
    ## phasing 5'
    trainingA = pd.merge(trainingDf, self.posOffsets, left_on=["chrom", "start"], right_on=["chrom", "pos"])
    trainingB = pd.merge(trainingDf, self.negOffsets, left_on=["chrom", "end"], right_on=["chrom", "pos"])
    trainingDf = pd.concat([trainingA, trainingB])
    trainingDf.rename(columns={'offset': 'five_offset'}, inplace=True)
    ## phasing 3'
    trainingA = pd.merge(trainingDf, self.posOffsets, left_on=["chrom", "end"], right_on=["chrom", "pos"])
    trainingB = pd.merge(trainingDf, self.negOffsets, left_on=["chrom", "start"], right_on=["chrom", "pos"])
    trainingDf = pd.concat([trainingA, trainingB])
    trainingDf.rename(columns={'offset': 'three_offset'}, inplace=True)
    ## filter a read by whether it has a-site that satisfies [9,18]
    if not self.RelE:
        trainingDf = trainingDf[((trainingDf['asite'] >= 9) & (trainingDf['asite'] <= 18))]
        trainingDf = trainingDf[(trainingDf['asite'] >= trainingDf['read_length'] / 2 - 1)]
    else:
        trainingDf = trainingDf[((trainingDf['asite'] >= 1) & (trainingDf['asite'] <= 8))]
    ## slice the dataframe to the variables needed for training data, removed "start_nt", "end_nt"
    trainingDf = trainingDf[["asite", "read_length", "five_offset", "three_offset", "gene_strand"]]
    trainingDf.to_csv(path_or_buf=self.outBam + '.training.txt', sep='\t', header=True, index=False)
    '''
    ## use dask to do merging
    #Construct a dask objects from a pandas objects
    trainingDf = dd.from_pandas(trainingDf, npartitions=2)
    trainingDf = dd.merge(trainingDf, self.offsets, left_on=["gene_name", "start"], right_on=["gene_name","pos"])
    trainingDf = trainingDf.rename(columns={'offset':'five_offset'})
    trainingDf = dd.merge(trainingDf, self.offsets, left_on=["gene_name", "end"], right_on=["gene_name","pos"])
    trainingDf = trainingDf.rename(columns={'offset':'three_offset'})
    trainingDf = trainingDf[((trainingDf['asite'] >= 12) & (trainingDf['asite'] <= 18))]
    trainingDf = trainingDf[(trainingDf['asite'] >= trainingDf['read_length'] / 2 - 1)]
    trainingDf = trainingDf[["asite", "read_length", "five_offset", "three_offset", "gene_strand"]].compute()
    trainingDf.to_csv(self.outBam + '.asite.txt', sep='\t', header=True, index=False)
    '''
def makeTestingData(self):
    ## create pandas df from bedtools intersect
    self.bedtool = pbt.BedTool(self.outBam + '.bed')
    testingData = self.bedtool.intersect(self.orf, wa=True, wb=True, sorted=True)
    time = str(datetime.now())
    print("[status]\tFinished intersecting the bedtool with ORFs:", time, flush=True)
    testingDf = testingData.to_dataframe(names=['chrom', 'start', 'end', 'name', 'read_length', 'strand', 'thickStart',
                                                'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts',
                                                'gene_chrom', 'gene_start', 'gene_end', 'gene_name', 'gene_score',
                                                'gene_strand', 'gene_thickStart', 'gene_thickEnd', 'gene_itemRgb',
                                                'gene_blockCount', 'gene_blockSizes', 'gene_blockStarts'],
                                         dtype={'blockSizes':'object','blockStarts':'object',
                                                'gene_blockSizes':'object','gene_blockStarts':'object'})
    ## phasing 5'
    testingA = pd.merge(testingDf, self.posOffsets, left_on=["chrom", "start"], right_on=["chrom","pos"])
    testingB = pd.merge(testingDf, self.negOffsets, left_on=["chrom", "end"], right_on=["chrom","pos"])
    testingDf = pd.concat([testingA, testingB])
    testingDf.rename(columns={'offset':'five_offset'}, inplace=True)
    ## phasing 3'
    testingA = pd.merge(testingDf, self.posOffsets, left_on=["chrom", "end"], right_on=["chrom","pos"])
    testingB = pd.merge(testingDf, self.negOffsets, left_on=["chrom", "start"], right_on=["chrom","pos"])
    testingDf = pd.concat([testingA, testingB])
    testingDf.rename(columns={'offset':'three_offset'}, inplace=True)
    ## slice the dataframe to the variables needed for training data
    testingDf = testingDf[["read_length", "five_offset", "three_offset", "gene_strand", "chrom", "start", "end"]]
    testingDf.to_csv(path_or_buf=self.outBam + '.testing.txt', sep='\t', header=True, index=False)
    '''
    ## use dask to do merging
    #Construct a dask objects from a pandas objects
    testingDf = dd.from_pandas(testingDf, npartitions=2)
    testingDf = dd.merge(testingDf, self.offsets, left_on=["gene_name", "start"], right_on=["gene_name","pos"])
    testingDf = testingDf.rename(columns={'offset':'five_offset'})
    testingDf = dd.merge(testingDf, self.offsets, left_on=["gene_name", "end"], right_on=["gene_name","pos"])
    testingDf = testingDf.rename(columns={'offset':'three_offset'})
    testingDf = testingDf[["read_length", "five_offset", "three_offset", "gene_strand", "chrom", "start", "end"]].compute()
    testingDf.to_csv(self.outBam + '.cds.txt', sep='\t', header=True, index=False)
    '''

#@profile
def genericGLM(self):
    y, X = dmatrices('ribosome_count ~ C(gene_name) + C(codon) + avg_prob', self.df)
    print(y)
    print("\n")
    print(X)

    #mod = NBin(y, X)
    #res = mod.fit('lbfgs')
    #print(res)


    ## alternative fit a gamma dist
    # formula = 'norm_count ~ C(gene_name) + C(codon)  + pair_prob '
    # print("[status]\tFormula: " + str(formula), flush=True)
    # mod = smf.glm(formula, self.df, family=sm.families.Gamma(link=sm.families.links.log))


class NBin(GenericLikelihoodModel):
    def __init__(self, endog, exog, **kwds):
        super(NBin, self).__init__(endog, exog, **kwds)

    def nloglikeobs(self, params):
        alph = params[-1]
        beta = params[:-1]
        ll = _ll_nb2(self.endog, self.exog, beta, alph)
        return -ll
    def fit(self, start_params=None, maxiter=1000, maxfun=5000, **kwds):
        if start_params == None:
            # Reasonable starting values
            start_params = np.append(np.zeros(self.exog.shape[1]), .5)
            start_params[0] = np.log(self.endog.mean())
        return super(NBin, self).fit(start_params=start_params,
                                     maxiter=maxiter, maxfun=maxfun,
                                     **kwds)

def rfFit(self):
    # column names
    self.colNames = list(self.traningDf.columns.values)
    self.colNames.remove("asite")
    self.X = np.array(pd.get_dummies(self.traningDf[self.colNames]))
    self.y = np.array(self.traningDf["asite"])
    ## feature selection
    self.clf = RandomForestClassifier(max_features=None, n_jobs=-1)
    self.clf = self.clf.fit(self.X, self.y)
    self.importances = self.clf.feature_importances_
    self.selector = RFECV(self.clf, step=1, cv=5)
    self.selector = self.selector.fit(self.X, self.y)
    self.sltX = self.selector.transform(self.X)
    print("[status]\tOptimal number of features based on recursive selection: %d" % self.selector.n_features_)
    ## define a new classifier for reduced features
    self.reducedClf = RandomForestClassifier(max_features=None, n_jobs=-1)
    self.reducedClf = self.reducedClf.fit(self.sltX, self.y)
    ## feature selection based on thresholding
    #cutoff = 0.02
    #self.selector = SelectFromModel(self.reducedClf, prefit=True, threshold=cutoff)
    #self.sltX = self.selector.transform(self.sltX)
    #numSltFeatures = self.sltX.shape[1]
    #print("[status]\t", str(numSltFeatures), "features with importance higher than", str(cutoff), flush=True)
    ## define a new classifier for reduced features
    #self.reducedClf = RandomForestClassifier(max_features=None, n_jobs=-1)
    #self.reducedClf = self.reducedClf.fit(self.sltX, self.y)
    ## cross validation
    scores = cross_val_score(self.reducedClf, self.sltX, self.y, cv=10)
    print("[status]\tAccuracy: %0.3f (+/- %0.3f)" % (scores.mean(), scores.std() * 2), flush=True)

def filterRegion(self):
    # create bedtool, filter and sort
    #self.startCodons = pbt.BedTool(self.startCodons).filter(lambda x: x.chrom != 'chrM')
    #self.startCodons = self.startCodons.sort()
    #self.orf = pbt.BedTool(self.orf).filter(lambda x: x.chrom != 'chrM')
    #self.orf = self.orf.sort()
    self.startCodons = pbt.BedTool(self.startCodons).sort()
    self.orf = pbt.BedTool(self.orf).sort()
    # find overlapping regions
    distinctStartCodons = self.startCodons.merge(c="1,2", o="count").filter(lambda x: int(x[4]) == 1)
    distinctOrfs = self.orf.merge(c="1,2", o="count").filter(lambda x: int(x[4]) == 1)
    distinctStartCodons = distinctStartCodons.intersect(distinctOrfs, wa=True, sorted=True)
    # filter start codon
    startCodonHash = set([(i[0], i[1], i[2]) for i in distinctStartCodons])
    self.startCodons = self.startCodons.filter(lambda x: (x[0], x[1], x[2]) in startCodonHash)


    def checkrpy(self):
        """

        :return:
        """
        # import R's utility package
        utils = rpackages.importr('utils')
        # select a mirror for R packages
        utils.chooseCRANmirror(ind=1)  # select the first mirror in the list
        # R package names
        packnames = ['glmnet']
        # R vector of strings
        from rpy2.robjects.vectors import StrVector
        # Selectively install what needs to be install.
        # We are fancy, just because we can.
        names_to_install = [x for x in packnames if not rpackages.isinstalled(x)]
        if len(names_to_install) > 0:
            utils.install_packages(StrVector(packnames), dependencies=True, repos='http://cran.rstudio.com/')

    def rpyFit(self):
        """

        :return:
        """
        n2r.activate()
        r = ro.r
        r.library('glmnet')
        #
        # from rpy2.robjects import pandas2ri
        # pandas2ri.activate()
        #
        self.df = pd.read_table("filtered.txt",  header=0)
        # offsets = ro.FloatVector(list(self.df.logTPM_scaled.transpose()))
        X = np.array(pd.get_dummies(self.df[["gene", "codon"]]).to_sparse(fill_value=0))
        y = np.array(self.df["ribosome_count"])
        y = ro.IntVector(list(y.transpose()))  # use factors
        fit = r['cv.glmnet'](X, y, family="poisson", alpha=0, parallel=True)
        # fit = r['cv.glmnet'](X, y, family="poisson", alpha=0, offset=offsets, parallel=True)
        fit = np.asanyarray(fit.rx2('lambda'))
        print(fit)
        # lambdas = [1600, 400, 200] + list(range(100, 0, -1))

