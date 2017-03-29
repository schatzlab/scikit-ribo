#!/usr/bin/env python

# ----------------------------------------
# scikit-ribo
# ----------------------------------------
# a module for a-site prediction
# ----------------------------------------
# author: Han Fang
# contact: hanfang.cshl@gmail.com
# website: hanfang.github.io
# date: 10/28/2016
# ----------------------------------------

from __future__ import print_function, division
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
import os
import sys
import argparse
import pandas as pd
import seaborn as sns
import numpy as np
import pybedtools as pbt
from sklearn import preprocessing, svm, tree
from sklearn.metrics import roc_curve, auc
from sklearn.preprocessing import label_binarize
from sklearn.multiclass import OneVsRestClassifier
from sklearn.model_selection import GridSearchCV, train_test_split, cross_val_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectFromModel, RFECV


class VisualizeAsite(object):
    ''' plotting a-site location distribution
    '''
    def __init__(self, training=None, RelE=None, output=None):
        self.training = training
        self.RelE = RelE
        self.output = output

    def plot(self):
        sns.set(font_scale=2)
        g0 = sns.FacetGrid(self.training, row="5_offset", col="read_length", margin_titles=True,
                           col_order= list(range(10,36)), row_order= list(range(3)))
        bins = np.linspace(12, 21, 10) if not self.RelE else np.linspace(1, 10, 10)
        g0.map(plt.hist, "asite", color="steelblue", bins=bins, lw=0,normed=True)
        if not self.RelE:
            g0.set(xticks=[12, 15, 18, 21])
        else:
            g0.set(xticks=[1, 4, 7, 10])
        plt.gcf()
        plt.savefig(self.output + "/" + "asite_5offset.pdf" )
        plt.clf()
        g1 = sns.FacetGrid(self.training, row="3_offset", col="read_length", margin_titles=True,
                           col_order= list(range(10,36)), row_order= list(range(3)))
        g1.map(plt.hist, "asite", color="steelblue", bins=bins, lw=0,normed=True)
        if not self.RelE:
            g1.set(xticks=[12, 15, 18, 21])
        else:
            g1.set(xticks=[1, 4, 7, 10])
        plt.gcf()
        plt.savefig(self.output + "/" + "asite_3offset.pdf" )
        plt.clf()


class PredictAsite(object):
    ''' model training - a-site prediction
    '''
    def __init__(self, training=None, cds=None, classifier="rf", RelE=None, pre=None, output=None, directory=None):
        self.training = training
        self.cds = cds
        self.cdsIdxFn = directory + "/" + pre + ".codons.df"
        self.classifier = classifier
        self.RelE = RelE
        self.output = output
        self.pre = pre

    def rfFit(self):
        # column names
        self.colNames = list(self.training.columns.values)
        self.colNames.remove("asite")
        self.X = np.array(pd.get_dummies(self.training[self.colNames]))
        self.y = np.array(self.training["asite"])
        ## feature selection
        self.clf = RandomForestClassifier(max_features=None, n_jobs=-1)
        self.clf = self.clf.fit(self.X, self.y)
        self.importances = self.clf.feature_importances_
        self.selector = RFECV(self.clf, step=1, cv=5)
        self.selector = self.selector.fit(self.X, self.y)
        self.sltX = self.selector.transform(self.X)
        print("[result]\tOptimal number of features by recursive selection: %d" % self.selector.n_features_, flush=True)
        ## define a new classifier for reduced features
        self.reducedClf = RandomForestClassifier(max_features=None, n_jobs=-1)
        self.reducedClf = self.reducedClf.fit(self.sltX, self.y)
        ## cross validation
        scores = cross_val_score(self.reducedClf, self.sltX, self.y, cv=10)
        print("[result]\tAccuracy: %0.3f (+/- %0.3f)" % (scores.mean(), scores.std() * 2), flush=True)

    def rfImportance(self):
        ## compute the std and index for the feature importance
        std = np.std([tree.feature_importances_ for tree in self.clf.estimators_], axis=0)
        idx = np.argsort(self.importances)[::-1]
        featureNames =(pd.get_dummies(self.training[self.colNames]).columns.values)
        importantFeatures = featureNames[idx]
        ## Plot the feature importances of the classifier
        plt.figure()
        plt.title("Feature importances")
        plt.bar(range(self.X.shape[1]), self.importances[idx], color=sns.xkcd_rgb["denim blue"], yerr=std[idx], align="center")
        plt.xticks(range(self.X.shape[1]), importantFeatures, rotation='vertical')
        plt.xlim([-1, 10])
        plt.ylim([0, 1])
        #plt.gca().tight_layout()
        plt.gcf()
        plt.savefig(self.output + "/" + "asite_feature_importances.pdf", facecolor="white")

    def rfPredict(self):
        ## create df for cds
        #self.cds = pd.read_table(self.cdsFn + ".txt",  header=0)
        cdsX = np.array(pd.get_dummies(self.cds[self.colNames]) )
        ## selected a subset of features and predict a-site
        sltcdsX = self.selector.transform(cdsX)
        self.cds["asite"] = self.reducedClf.predict(sltcdsX)

    def svmFit(self):
        ## grid search
        self.clf = svm.SVC()
        paramGrid = [{'C': [ 0.01, 0.1, 1, 10, 100, 1000, 10000]}]
        self.clfGs = GridSearchCV(estimator=self.clf, param_grid=paramGrid, n_jobs=-1)
        self.clfGs.fit(self.X, self.y)
        print("[result]\t best estimator parameters: c=", self.clfGs.best_estimator_.C, flush=True)
        ## model fitting and cross validation
        self.clf = svm.SVC(C=self.clfGs.best_estimator_.C)
        scores = cross_val_score(self.clf, self.X, self.y, cv=10)
        print("[result]\tAccuracy: %0.3f (+/- %0.3f)" % (scores.mean(), scores.std() * 2), flush=True)

    def rocCurve(self):
        ''' plotting multi-class roc curve
        '''
        # shuffle and split training and test sets
        clf = self.reducedClf if self.classifier == "rf" else self.clf
        self.OvrClf = OneVsRestClassifier(clf)
        classes = list(range(9,19)) if not self.RelE else list(range(1,9))
        self.y = label_binarize(self.y, classes=classes)
        nClasses = self.y.shape[1]
        X_train, X_test, y_train, y_test = train_test_split(self.X, self.y, test_size=.5, random_state=0)
        if self.classifier == "rf":
            y_score = self.OvrClf.fit(X_train, y_train).predict_proba(X_test)
        else:
            y_score = self.OvrClf.fit(X_train, y_train).decision_function(X_test)
        # Compute ROC curve and ROC area for each class
        fpr, tpr, roc_auc = {}, {}, {}
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
            pos = classes[i]
            plt.plot(fpr[i], tpr[i], label='A-site @ {0} (area = {1:0.2f})'
                     ''.format(pos, roc_auc[i]))
        #
        plt.plot([0, 1], [0, 1], 'k--')
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate', fontsize=18)
        plt.ylabel('True Positive Rate', fontsize=18)
        plt.tick_params(axis='both', which='major', labelsize=18)
        plt.legend(loc="lower right",fontsize=12)
        plt.gcf()
        plt.savefig(self.output + "/" + "asite_roc.pdf")

    def recoverAsite(self):
        ## adjust by the a-site location and calculate the a-site location in nt space, -1 is the missing value
        if not self.RelE:
            self.cds['a_start'] = np.where(self.cds['gene_strand'] == '+', (self.cds['start'] + self.cds['asite']),
                                           (-1)).astype(int)
            self.cds['a_end'] = np.where(self.cds['gene_strand'] == '+', (self.cds['a_start'] + 3),
                                         (self.cds['end'] - self.cds['asite'])).astype(int)
            self.cds['a_start'] = np.where(self.cds['gene_strand'] == '-', (self.cds['a_end'] - 3),
                                           (self.cds['a_start'])).astype(int)
        else:
            self.cds['a_start'] = np.where(self.cds['gene_strand'] == '+', (self.cds['end'] - self.cds['asite']),
                                           (-1)).astype(int)
            self.cds['a_end'] = np.where(self.cds['gene_strand'] == '+', (self.cds['a_start'] + 3),
                                         (self.cds['start'] + self.cds['asite'])).astype(int)
            self.cds['a_start'] = np.where(self.cds['gene_strand'] == '-', (self.cds['a_end'] - 3),
                                           (self.cds['a_start'])).astype(int)
        # remove start/end for reads
        self.cds.drop(['start', 'end'], axis=1, inplace=True)
        ## use to group by command to retrieve ribosome coverage
        cnt = self.cds.groupby(["chrom", "a_start", "a_end", "strand"])
        cnt = cnt.size().reset_index(name="ribosome_count")
        ## left outer join the null df and the groupby_df_count to get ribsome counts at each position
        cdsIdx = pd.read_table(self.cdsIdxFn, header=0)
        riboCnt = pd.merge(cdsIdx, cnt, how="left", left_on=["chrom", "start", "end", "gene_strand"],
                           right_on=["chrom", "a_start", "a_end", "strand"])
        riboCnt.drop(['a_start', 'a_end', 'strand'], axis=1, inplace=True)
        riboCnt["ribosome_count"].fillna(value=0, inplace=True)
        riboCnt["ribosome_count"] = riboCnt["ribosome_count"].astype(int)
        riboCnt = riboCnt.sort_values(by=["chrom", "start", "end"])
        return riboCnt


# ----------------------------------------
#  main
# ----------------------------------------
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help="input folder, required")
    parser.add_argument("-p", help="prefix of index for the CDS region, required")
    parser.add_argument("-c", help="classifier to use, random forest (rf) or svm, optional, default: rf", default="rf")
    parser.add_argument("-e", help="whether the sample involved RelE, Default: F", default='F', type=str)
    # check if there is any argument
    if len(sys.argv) <= 1:
        parser.print_usage()
        sys.exit(1)
    else:
        args = parser.parse_args()
    # process the file if the input files exist
    if (args.i != None and args.p != None):
        sys.stderr.write("[status]\tprocessing the input file: " + args.i + "\n")
        input = args.i
        asite_fn = input + "/riboseq.training"
        cds_fn = input + "/riboseq.cds"
        cds_idx = input + "/" + args.p + ".codons.df"
        classifier = args.c
        RelE = False if args.e == 'F' else True
        sys.stderr.write("[execute]\tplotting the a-site location distribution from " + str(asite_fn) + "\n")
        asite = visualize_asite(asite_fn, RelE)
        asite.plot()
        sys.stderr.write("[execute]\tstart the process of a-site prediction" + "\n")
        # model training
        model = predict_site(asite_fn, cds_fn, cds_idx, classifier, RelE)
        if classifier == "rf":
            sys.stderr.write("[execute]\tperform model training and cross validation on the training data" + "\n")
            model.rfFit()
            sys.stderr.write("[execute]\tplotting the bar plot of the feature importance" + "\n")
            model.rfImportance()
            sys.stderr.write("[execute]\tplot roc curve based on cross validation" + "\n")
            model.rocCurve()
            sys.stderr.write("[execute]\tpredicting the a-site from the cds regions" + "\n")
            model.rfPredict()
            sys.stderr.write("[execute]\tlocalize the a-site codon and create coverage df" + "\n")
            model.recoverAsite()
        elif classifier == "svm":
            sys.stderr.write("[execute]\tperform SVM classifier training" + "\n")
            model.svmFit()
            sys.stderr.write("[execute]\tplot roc curve based on cross validation" + "\n")
            model.rocCurve()
        # end
        sys.stderr.write("[status]\tA-site module finished" + "\n")
    else:
        sys.stderr.write("[error]\tmissing argument" + "\n")
        parser.print_usage()
