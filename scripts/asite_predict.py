#!/usr/bin/env python

## ----------------------------------------
## scikit-ribo
## ----------------------------------------
## a module for preprocessing bam files
## ----------------------------------------
## author: Han Fang
## contact: hanfang.cshl@gmail.com
## website: hanfang.github.io
## date: 1/28/2016
## ----------------------------------------

from __future__ import print_function, division
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import sys
import pandas as pd
import seaborn as sns
import numpy as np
import argparse
from sklearn import linear_model, cross_validation, preprocessing, svm, tree
from sklearn.metrics import roc_curve, auc
from sklearn.cross_validation import train_test_split
from sklearn.preprocessing import label_binarize
from sklearn.multiclass import OneVsRestClassifier
from sklearn.grid_search import GridSearchCV

class VisualizeAsite:
    ''' define a class for plotting a-site location distribution
    '''
    def __init__(self, asite_fn):
        self.asite_fn = asite_fn
        self.asite_df = pd.read_table(self.asite_fn,  header=0)

    def plot(self):
        sns.set(font_scale=2)
        g0 = sns.FacetGrid(self.asite_df,
                           row="offset",
                           col="read_length",
                           margin_titles=True,
                           col_order= list(range(25,36)),
                           row_order= list(range(3)))
        bins = np.linspace(12, 21, 10)
        g0.map(plt.hist, "a_site", color="steelblue", bins=bins, lw=0,normed=True)
        g0.set(xticks=[12,15,18,21])
        plt.gcf()
        plt.savefig( asite_fn +".asite.pdf" )
        plt.clf()

class TrainModel:
    ''' define a class for model training - a-site prediction
    '''
    def __init__(self, object):
        self.asite_df = object.asite_df
        self.df_colnames = list(object.asite_df.columns.values)
        self.df_colnames.remove("a_site")
        self.X = np.array( pd.get_dummies( object.asite_df[self.df_colnames] ))
        self.y = np.array(object.asite_df["a_site"])

    def fit(self):
        poly = preprocessing.PolynomialFeatures(2)
        poly.fit_transform(self.X)
        # clf = svm.SVC()
        
        ## grid search
        self.clf = linear_model.LogisticRegression()
        param_grid = [ { 'C': [ 0.01, 0.1, 1, 10, 100, 1000, 10000] } ]
        self.clf_gs = GridSearchCV(estimator=self.clf, param_grid=param_grid, n_jobs=-1)
        self.clf_gs.fit(self.X, self.y)
        print("[status]\t best estimator parameters: c=", self.clf_gs.best_estimator_.C)

        ## model fitting and cross validation
        self.clf = linear_model.LogisticRegression( C=self.clf_gs.best_estimator_.C,
                                                    solver="lbfgs"
                                                    )
        
        # self.clf = linear_model.LogisticRegressionCV(n_jobs=-1, cv=10)
        scores = cross_validation.cross_val_score(self.clf, self.X, self.y, cv=10)
        print( self.clf.get_params() )
        print("Accuracy: %0.3f (+/- %0.3f)" % (scores.mean(), scores.std() * 2))

    def roc_curve(self):
        ''' define a class for plotting multi-class roc curve
        '''

        # shuffle and split training and test sets
        self.OVR_clf = OneVsRestClassifier(self.clf)
        self.y = label_binarize(self.y, classes=[12,13,14,15,16,17,18])
        n_classes = self.y.shape[1]
        X_train, X_test, y_train, y_test = train_test_split(self.X, self.y, test_size=.5, random_state=0)
        y_score = self.OVR_clf.fit(X_train, y_train).decision_function(X_test)

        # Compute ROC curve and ROC area for each class
        fpr = dict()
        tpr = dict()
        roc_auc = dict()
        for i in range(n_classes):
            fpr[i], tpr[i], _ = roc_curve(y_test[:, i], y_score[:, i])
            roc_auc[i] = auc(fpr[i], tpr[i])

        # Compute micro-average ROC curve and ROC area
        fpr["micro"], tpr["micro"], _ = roc_curve(y_test.ravel(), y_score.ravel())
        roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])

        # Plot ROC curve
        sns.reset_orig()
        plt.clf()
        plt.figure()
        plt.plot(fpr["micro"], tpr["micro"],'--', linewidth=3,
                 label='micro-average (area = {0:0.2f})'
                 ''.format(roc_auc["micro"]))
        for i in range(n_classes):
            plt.plot(fpr[i], tpr[i], label='A-site @ {0} (area = {1:0.2f})'
                     ''.format(i+12, roc_auc[i]))

        plt.plot([0, 1], [0, 1], 'k--')
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate', fontsize=18)
        plt.ylabel('True Positive Rate', fontsize=18)
        plt.tick_params(axis='both', which='major', labelsize=18)
        plt.legend(loc="lower right",fontsize=12)
        plt.gcf()
        plt.savefig(asite_fn + ".roc.pdf")

## ----------------------------------------
## the main work
## ----------------------------------------
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help="input a-site file")
    #parser.add_argument("-q", help="minimum mapq allowed", default=10, type=int)

    ## check if there is any argument
    if len(sys.argv) <= 1:
        parser.print_usage()
        sys.exit(1)
    else:
        args = parser.parse_args()

    ## process the file if the input files exist
    if (args.i != None):
        print("[status]\tprocessing the input file: " + args.i)
        asite_fn = args.i

        print("[execute]\tplotting the a-site location distribution from " + str(asite_fn))
        asite_loc = VisualizeAsite(asite_fn)
        asite_loc.plot()

        print("[execute]\tperform cross validation on the training data")
        model = TrainModel(asite_loc)
        model.fit()
        model.roc_curve()

    else:
        print("[error]\tmissing argument")
        parser.print_usage()
