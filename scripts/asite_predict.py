#!/usr/bin/env python

## ----------------------------------------
## scikit-ribo
## ----------------------------------------
## a module for a-site prediction
## ----------------------------------------
## author: Han Fang
## contact: hanfang.cshl@gmail.com
## website: hanfang.github.io
## date: 10/28/2016
## ----------------------------------------

from __future__ import print_function, division
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
import os
import sys
import pandas as pd
import seaborn as sns
import numpy as np
import pybedtools as pbt
import argparse
from sklearn import preprocessing, svm, tree
from sklearn.metrics import roc_curve, auc
from sklearn.preprocessing import label_binarize
from sklearn.multiclass import OneVsRestClassifier
from sklearn.model_selection import GridSearchCV, train_test_split, cross_val_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectFromModel, RFECV


class visualizeAsite(object):
    ''' define a class for plotting a-site location distribution
    '''
    def __init__(self, asiteFn):
        self.asiteFn = asiteFn
        self.asiteDf = pd.read_table(self.asiteFn,  header=0)

    def plot(self):
        sns.set(font_scale=2)
        g0 = sns.FacetGrid(self.asiteDf, row="five_offset", col="read_length", margin_titles=True,
                           col_order= list(range(25,36)),
                           row_order= list(range(3)))
        bins = np.linspace(12, 21, 10)
        g0.map(plt.hist, "asite", color="steelblue", bins=bins, lw=0,normed=True)
        g0.set(xticks=[12,15,18,21])
        plt.gcf()
        plt.savefig(self.asiteFn +".asite_five_offset.pdf" )
        plt.clf()

        g1 = sns.FacetGrid(self.asiteDf, row="three_offset", col="read_length", margin_titles=True,
                           col_order= list(range(25,36)),
                           row_order= list(range(3)))
        bins = np.linspace(12, 21, 10)
        g1.map(plt.hist, "asite", color="steelblue", bins=bins, lw=0,normed=True)
        g1.set(xticks=[12,15,18,21])
        plt.gcf()
        plt.savefig(self.asiteFn +".asite_three_offset.pdf" )
        plt.clf()

class trainModel(object):
    ''' define a class for model training - a-site prediction
    '''
    def __init__(self, object, testingFn, cdsIdxFn, classifier):
        self.traningDf = object.asiteDf
        self.asiteFn = object.asiteFn
        self.colNames = list(self.traningDf.columns.values)
        self.colNames.remove("asite")
        self.X = np.array(pd.get_dummies(self.traningDf[self.colNames]) )
        self.y = np.array(self.traningDf["asite"])
        self.testingFn = testingFn
        self.testingDf = pd.read_table(self.testingFn, header=0) # .dropna()
        self.cdsIdxFn = cdsIdxFn
        self.cdsIdxDf = pd.read_table(self.cdsIdxFn, header=0)
        self.featureNames =(pd.get_dummies(self.traningDf[self.colNames]).columns.values)
        self.classifier = classifier

    def rfFit(self):
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
        cutoff = 0.01
        self.selector = SelectFromModel(self.reducedClf, prefit=True, threshold=cutoff)
        self.sltX = self.selector.transform(self.sltX)
        numSltFeatures = self.sltX.shape[1]
        print("[status]\t", str(numSltFeatures), "features with importance higher than", str(cutoff), flush=True)
        ## define a new classifier for reduced features
        self.reducedClf = RandomForestClassifier(max_features=None, n_jobs=-1)
        self.reducedClf = self.reducedClf.fit(self.sltX, self.y)
        ## cross validation
        scores = cross_val_score(self.reducedClf, self.sltX, self.y, cv=10)
        print("[status]\tAccuracy: %0.3f (+/- %0.3f)" % (scores.mean(), scores.std() * 2), flush=True)

    def rfImportance(self):
        ## compute the std and index for the feature importance
        std = np.std([tree.feature_importances_ for tree in self.clf.estimators_],    axis=0)
        idx = np.argsort(self.importances)[::-1]
        importantFeatures = self.featureNames[idx]

        ## Plot the feature importances of the classifier
        plt.figure()
        plt.title("Feature importances")
        plt.bar(range(self.X.shape[1]), self.importances[idx], color=sns.xkcd_rgb["denim blue"], yerr=std[idx], align="center")
        plt.xticks(range(self.X.shape[1]), importantFeatures, rotation='vertical')
        plt.xlim([-1, 10])
        plt.ylim([0, 1])
        #plt.gca().tight_layout()
        plt.gcf()
        plt.savefig(self.asiteFn + ".feature_importances.pdf", facecolor="white")

    def rfPredict(self):
        ## create df for testing data
        self.testingDf = pd.read_table(self.testingFn,  header=0)
        #colNames = list(testingDf.columns.values)
        #namesToExclude = set(["gene_start", "gene_end", "gene_name", "chrom", "start", "end", "name", "mapq", "strand"])
        #remainColNames = [x for x in colNames if x not in namesToExclude]
        #testingX = np.array(pd.get_dummies(testingDf[remainColNames]) )
        testingX = np.array(pd.get_dummies(self.testingDf[self.colNames]) )
        ## selected a subset of features and predict a-site
        sltTestingX = self.selector.transform(testingX)
        self.testingDf["asite"] = self.reducedClf.predict(sltTestingX)
        #self.testingDfOut = testingDf[["chrom", "start", "end", "name", "strand", "asite", "gene_name", "gene_strand"]
        #                              + self.colNames ]
        self.testingDf.to_csv(path_or_buf=self.testingFn + '.predicted.txt', sep='\t', header=True, index=False)

    def svmFit(self):
        ## grid search
        self.clf = svm.SVC()
        paramGrid = [{'C': [ 0.01, 0.1, 1, 10, 100, 1000, 10000]}]
        self.clfGs = GridSearchCV(estimator=self.clf, param_grid=paramGrid, n_jobs=-1)
        self.clfGs.fit(self.X, self.y)
        print("[status]\t best estimator parameters: c=", self.clfGs.best_estimator_.C, flush=True)
        ## model fitting and cross validation
        self.clf = svm.SVC(C=self.clfGs.best_estimator_.C)
        scores = cross_val_score(self.clf, self.X, self.y, cv=10)
        print("[status]\tAccuracy: %0.3f (+/- %0.3f)" % (scores.mean(), scores.std() * 2), flush=True)

    def rocCurve(self):
        ''' define a class for plotting multi-class roc curve
        '''
        # shuffle and split training and test sets
        clf = self.reducedClf if self.classifier == "rf" else self.clf
        self.OvrClf = OneVsRestClassifier(clf)
        self.y = label_binarize(self.y, classes=[12,13,14,15,16,17,18])
        nClasses = self.y.shape[1]
        X_train, X_test, y_train, y_test = train_test_split(self.X, self.y, test_size=.5, random_state=0)
        if self.classifier == "rf":
            y_score = self.OvrClf.fit(X_train, y_train).predict_proba(X_test)
        else:
            y_score = self.OvrClf.fit(X_train, y_train).decision_function(X_test)
        # Compute ROC curve and ROC area for each class
        fpr = dict()
        tpr = dict()
        roc_auc = dict()
        for i in range(nClasses):
            fpr[i], tpr[i], _ = roc_curve(y_test[:, i], y_score[:, i])
            roc_auc[i] = auc(fpr[i], tpr[i])
        # Compute micro-average ROC curve and ROC area
        fpr["micro"], tpr["micro"], _ = roc_curve(y_test.ravel(), y_score.ravel())
        roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])
        # Plot ROC curve
        sns.reset_orig()
        plt.clf()
        plt.figure()
        plt.plot(fpr["micro"], tpr["micro"],'--', linewidth=3, label='micro-average (area = {0:0.2f})'
                 ''.format(roc_auc["micro"]))
        for i in range(nClasses):
            plt.plot(fpr[i], tpr[i], label='A-site @ {0} (area = {1:0.2f})'
                     ''.format(i+12, roc_auc[i]))
        #
        plt.plot([0, 1], [0, 1], 'k--')
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate', fontsize=18)
        plt.ylabel('True Positive Rate', fontsize=18)
        plt.tick_params(axis='both', which='major', labelsize=18)
        plt.legend(loc="lower right",fontsize=12)
        plt.gcf()
        plt.savefig(self.asiteFn + ".roc.pdf")

    def recoverAsite(self):
        ## [temp] import the predicted a-site location df
        # self.testingDf = pd.read_table(self.testingFn+'.predicted.txt',  header=0)

        ## adjust by the a-site location and calculate the a-site location in nt space, -1 is the missing value
        self.testingDf['asite_start'] = np.where(self.testingDf['gene_strand'] == '+',
                                                 (self.testingDf['start'] + self.testingDf['asite']),
                                                 (-1) ).astype(int)
        self.testingDf['asite_end'] = np.where(self.testingDf['gene_strand'] == '+',
                                               (self.testingDf['asite_start'] + 3),
                                               (self.testingDf['end'] - self.testingDf['asite']) ).astype(int)
        self.testingDf['asite_start'] = np.where(self.testingDf['gene_strand'] == '-',
                                                 (self.testingDf['asite_end'] - 3),
                                                 (self.testingDf['asite_start']) ).astype(int)

        ## use to group by command to retrieve ribosome coverage
        coverageDf = self.testingDf.groupby(["chrom", "asite_start", "asite_end"])
        coverageDf = coverageDf.size().reset_index(name="ribosome_count")
        # groupbyDfCount.to_csv(path_or_buf='groupby_df_count.txt', sep='\t', header=True, index=False)

        ## left outer join the null df and the groupby_df_count to get ribsome counts at each position
        self.riboCountDf = pd.merge(self.cdsIdxDf, coverageDf, how = "left",
                                    on = ["chrom", "asite_start", "asite_end"])
        self.riboCountDf["ribosome_count"].fillna(value=0, inplace=True)
        self.riboCountDf["ribosome_count"] = self.riboCountDf["ribosome_count"].astype(int)
        self.riboCountDf = self.riboCountDf.sort_values(by=["chrom", "asite_start", "asite_end"])
        self.riboCountDf.to_csv(path_or_buf=self.testingFn + '.ribocount.txt', sep='\t', header=True, index=False)


## ----------------------------------------
## the main work
## ----------------------------------------
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help="input a-site file, required")
    parser.add_argument("-t", help="input testing data for the entire CDS region, required")
    parser.add_argument("-d", help="input index data-frame for the entire CDS region, required")
    parser.add_argument("-c", help="classifier to use, random forest (rf) or svm, optional, default: rf", default="rf")
    parser.add_argument("-o", help="output path, required")

    ## check if there is any argument
    if len(sys.argv) <= 1:
        parser.print_usage()
        sys.exit(1)
    else:
        args = parser.parse_args()

    ## process the file if the input files exist
    if (args.i != None and args.t != None):
        print("[status]\tprocessing the input file: " + args.i, flush=True)
        asite_fn = args.i
        cds_fn = args.t
        cds_idx = args.d
        classifier = args.c
        output = args.o

        cmd = "mkdir -p " + output
        os.system(cmd)

        print("[execute]\tplotting the a-site location distribution from " + str(asite_fn), flush=True)
        asite_loc = visualizeAsite(asite_fn)
        asite_loc.plot()

        print("[execute]\tstart the process of a-site prediction", flush=True)
        model = trainModel(asite_loc, cds_fn, cds_idx, classifier)

        if classifier == "rf":
            print("[execute]\tperform model training and cross validation on the training data", flush=True)
            model.rfFit()
            print("[execute]\tplotting the bar plot of the feature importance", flush=True)
            model.rfImportance()
            print("[execute]\tplot roc curve based on cross validation", flush=True)
            model.rocCurve()
            print("[execute]\tpredicting the a-site from the testing data", flush=True)
            model.rfPredict()
            print("[execute]\tlocalize the a-site codon and create coverage df", flush=True)
            model.recoverAsite()

        elif classifier == "svm":
            print("[execute]\tperform SVM classifier training", flush=True)
            model.svmFit()
            print("[execute]\tplot roc curve based on cross validation", flush=True)
            model.rocCurve()

        ## end
        print("[status]\tA-site module finished.", flush=True)

    else:
        print("[error]\tmissing argument", flush=True)
        parser.print_usage()
