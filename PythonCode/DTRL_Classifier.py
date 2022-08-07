import os
import sys
from setuptools.command.test import test
from sklearn import svm
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn import tree
from sklearn.ensemble import RandomForestClassifier, BaggingClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neural_network import MLPClassifier
from sklearn import preprocessing
from sklearn.decomposition import PCA, KernelPCA, IncrementalPCA
from joblib import dump, load
from os import path
import random
#import ROC
from sklearn import metrics



def printHeadingByDatasetNumber(i, outputFileName):
    outputFile = open(outputFileName, "a")
    if i == 1:
        print("\n==== Low Additive Dataset : ====", file=outputFile)
    elif i == 2:
        print("\n==== Low Mixed Dataset : =====", file=outputFile)
    elif i == 3:
        print("\n==== Low Replacing Dataset : ====", file=outputFile)
    elif i == 4:
        print("\n==== Medium Additive Dataset : ====", file=outputFile)
    elif i == 5:
        print("\n==== Medium Mixed Dataset: ====", file=outputFile)
    elif i == 6:
        print("\n==== Medium Replacing Dataset: ====", file=outputFile)
    elif i == 7:
        print("\n==== High Additive Dataset: ====", file=outputFile)
    elif i == 8:
        print("\n==== High Mixed Dataset: ==== ", file=outputFile)
    elif i == 9:
        print("\n==== High Replacing Dataset: ==== ", file=outputFile)


def printTrueLabelSuppliedResults(output, test_Y, i1, auc, outputFileName):
    printHeadingByDatasetNumber(i1, outputFileName)

    true_additive = 0.0
    misclassified_additive = 0.0
    false_additive = 0.0
    true_replacing = 0.0
    misclassified_replacing = 0.0
    false_replacing = 0.0

    for i in range(len(output)):
        if output[i] == 1 and test_Y[i] == 1:
            true_additive += 1
        elif test_Y[i] == 1:
            misclassified_additive += 1
        if output[i] == 1 and test_Y[i] != 1:
            false_additive += 1
        if output[i] == 0 and test_Y[i] == 0:
            true_replacing += 1
        elif test_Y[i] == 0:
            misclassified_replacing += 1
        if output[i] == 0 and test_Y[i] != 0:
            false_replacing += 1

    outputFile = open(outputFileName, "a")
    print("True additive: " + str(true_additive), file=outputFile)
    print("False negative for additive: " + str(misclassified_additive), file=outputFile)
    print("False positive for additive: " + str(false_additive) + "\n", file=outputFile)
    print("True replacing: " + str(true_replacing), file=outputFile)
    print("False negative for replacing: " + str(misclassified_replacing), file=outputFile)
    print("False positive for replacing: " + str(false_replacing), file=outputFile)

    additive_accuracy = true_additive / (true_additive + misclassified_additive + 0.00001) * 100.0
    replacing_accuracy = true_replacing / (true_replacing + misclassified_replacing + .00001) * 100.0
    false_positive_rate_additive = false_additive / (true_additive + false_additive + .00001) * 100.0
    false_positive_rate_replacing = false_replacing / (true_replacing + false_replacing + .00001) * 100.0
    false_negative_rate_additive = misclassified_additive / (true_additive + misclassified_additive + 0.00001) * 100.0
    false_negative_rate_replacing = misclassified_replacing / (
            true_replacing + misclassified_replacing + .00001) * 100.0
    if true_replacing < 5:
        replacing_accuracy = 0
        false_positive_rate_replacing = 0
        false_negative_rate_replacing = 0
    if true_additive < 5:
        additive_accuracy = 0
        false_positive_rate_additive = 0
        false_negative_rate_additive = 0

    tp = true_additive
    fp = misclassified_replacing
    tn = true_replacing
    fn = misclassified_additive
    precision = tp / (tp + fp + 0.0000001)
    recall = tp / (tp + fn + 0.0000001)
    f_score = 2 * (precision * recall) / (precision + recall + .0000001)
    print("\nAdditive accuracy: " + str(additive_accuracy) + " percent", file=outputFile)
    print("Replacing accuracy: " + str(replacing_accuracy) + " percent\n", file=outputFile)
    print("False positive rate for Additive transfers: " + str(false_positive_rate_additive) + " percent",
          file=outputFile)
    print("False negative rate for Additive transfers: " + str(false_negative_rate_additive) + " percent",
          file=outputFile)

    print("\nFalse positive rate for Replacing transfers: " + str(false_positive_rate_replacing) + " percent",
          file=outputFile)
    print("\nFalse negative rate for Replacing transfers: " + str(false_negative_rate_replacing) + " percent",
          file=outputFile)

    print("\nPrecision: " + str(precision), file=outputFile)
    print("\nRecall: " + str(recall), file=outputFile)
    print("\nF-measure: " + str(f_score), file=outputFile)
    print("\nAUC: " + str(auc), file=outputFile)
    print('\n')


def printRangerDependentResults(output, rangerLabel, test_Y, i, auc, outputFileName):
    printHeadingByDatasetNumber(i, outputFileName)

    true_additive = 0.0
    misclassified_additive = 0.0
    false_additive = 0.0
    true_replacing = 0.0
    misclassified_replacing = 0.0
    false_replacing = 0.0

    for i in range(len(output)):
        if output[i] == 1 and rangerLabel[i] != 2 and test_Y[i] == 1:
            true_additive += 1
        elif test_Y[i] == 1:
            misclassified_additive += 1
        if output[i] == 1 and rangerLabel[i] != 2 and test_Y[i] != 1:
            false_additive += 1
        if output[i] == 0 and rangerLabel[i] != 2 and test_Y[i] == 0:
            true_replacing += 1
        elif test_Y[i] == 0:
            misclassified_replacing += 1
        if output[i] == 0 and rangerLabel[i] != 2 and test_Y[i] != 0:
            false_replacing += 1

    outputFile = open(outputFileName, "a")
    print("True additive: " + str(true_additive), file=outputFile)
    print("False negative for additive: " + str(misclassified_additive), file=outputFile)
    print("False positive for additive: " + str(false_additive) + "\n", file=outputFile)
    print("True replacing: " + str(true_replacing), file=outputFile)
    print("False negative for replacing: " + str(misclassified_replacing), file=outputFile)
    print("False positive for replacing: " + str(false_replacing), file=outputFile)

    additive_accuracy = true_additive / (true_additive + misclassified_additive + 0.00001) * 100.0
    replacing_accuracy = true_replacing / (true_replacing + misclassified_replacing + .00001) * 100.0
    false_positive_rate_additive = false_additive / (true_additive + false_additive + .00001) * 100.0
    false_positive_rate_replacing = false_replacing / (true_replacing + false_replacing + .00001) * 100.0
    false_negative_rate_additive = misclassified_additive / (true_additive + misclassified_additive + 0.00001) * 100.0
    false_negative_rate_replacing = misclassified_replacing / (
            true_replacing + misclassified_replacing + .00001) * 100.0
    if true_replacing < 5:
        replacing_accuracy = 0
        false_positive_rate_replacing = 0
        false_negative_rate_replacing = 0
    if true_additive < 5:
        additive_accuracy = 0
        false_positive_rate_additive = 0
        false_negative_rate_additive = 0

    print("\nAdditive accuracy: " + str(additive_accuracy) + " percent", file=outputFile)
    print("Replacing accuracy: " + str(replacing_accuracy) + " percent\n", file=outputFile)
    print("False positive rate for Additive transfers: " + str(false_positive_rate_additive) + " percent",
          file=outputFile)
    print("False negative rate for Additive transfers: " + str(false_negative_rate_additive) + " percent",
          file=outputFile)

    print("\nFalse positive rate for Replacing transfers: " + str(false_positive_rate_replacing) + " percent",
          file=outputFile)
    print("\nFalse negative rate for Replacing transfers: " + str(false_negative_rate_replacing) + " percent",
          file=outputFile)
    print("\nAUC: " + str(auc), file=outputFile)
    print('\n')


def handle_ranger_dependant_ROC_curves(test_Y, rangerLabel, true_label_supplied_predictions):
    true_labels = []
    selected_predictions = []
    for i in range(len(test_Y)):
        if (rangerLabel[i] != 2 and test_Y[i] != 2):
            true_labels.append(test_Y[i])
            selected_predictions.append(true_label_supplied_predictions[i])
    ROC.plot_roc(true_labels, selected_predictions)


def classify(classifier, test_data, i, normalizer):
    test_X = test_data.loc[:, ['height',                             
                              'mappingCountHeuristic',                             
                              'lostGene2Heuristic',
                              'lostGene1Heuristic',
                              'lostGene3Heuristic',                             
                              'geneFrequencyHeuristic'
                              ]]
    test_Y = test_data.trueLabel
    lossHeuristicLabel = test_data.loc[:, "lossHeuristic"]
    currentHeuristicLabel = test_data.loc[:, "currentHeuristicLabel"]

    rangerLabel = test_data.loc[:, "currentHeuristicLabel"]

    #test_X[['number_of_descendents']] = scalar.fit_transform(test_X[['number_of_descendents']])
    test_X[['height']] = scalar.fit_transform(test_X[['height']])
    #print(test_X[['height']])
    output = classifier.predict(test_X)
    true_label_supplied_predictions = classifier.predict_proba(test_X)[:, 1]
    # ROC.plot_roc(test_Y, true_label_supplied_predictions)
    # handle_ranger_dependant_ROC_curves(test_Y, rangerLabel, true_label_supplied_predictions)

    fpr, tpr, thresholds = metrics.roc_curve(test_Y, true_label_supplied_predictions, pos_label=1)
    auc = metrics.auc(fpr, tpr)

    printTrueLabelSuppliedResults(output, test_Y, i, auc, "ML_output.txt")
    #printTrueLabelSuppliedResults(currentHeuristicLabel, test_Y, i,auc, "current_heuristic_output.txt")
    #printTrueLabelSuppliedResults(lossHeuristicLabel, test_Y, i,auc, "loss_heuristic_output.txt")
    printRangerDependentResults(output, rangerLabel, test_Y, i,auc, "ML_ranger_dependant_output.txt")


def get_classifier(training_data):
    X = training_data.loc[:, ['height',                             
                              'mappingCountHeuristic',                             
                              'lostGene2Heuristic',
                              'lostGene1Heuristic',
                              'lostGene3Heuristic',                             
                              'geneFrequencyHeuristic'
                              ]]
    Y = training_data.trueLabel

    # X = SelectKBest(chi2, k=4).fit_transform(X, Y)

    transformer = preprocessing.MinMaxScaler()

    X[['height']] = transformer.fit_transform(X[['height']])
    #X[['number_of_descendents']] = transformer.fit_transform(X[['number_of_descendents']])
    # pca.fit(X)
    # X = pca.transform(X)

    #  X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size = 0.2, random_state = 42)

    # clf = BaggingClassifier(KNeighborsClassifier(n_neighbors=3), max_samples=.8, max_features= .8)
    # clf = KNeighborsClassifier(n_neighbors=3)
    clf = RandomForestClassifier( n_estimators=10)
    # clf = svm.SVC()

    # clf = tree.ExtraTreeClassifier()
    ''' 
    clf = MLPClassifier(activation='relu', alpha=1e-05, batch_size='auto',
              beta_1=0.9, beta_2=0.999, early_stopping=False,
              epsilon=1e-08, hidden_layer_sizes=(10,5,2),
              learning_rate='constant', learning_rate_init=0.01,
             max_iter=1000, momentum=0.9, n_iter_no_change=10,
             nesterovs_momentum=True, power_t=0.5,  random_state=1,
             shuffle=True, solver='lbfgs', tol=0.0001,
              validation_fraction=0.1, verbose=False, warm_start=False)
    '''
    clf.fit(X, Y)

    # classify_train_test_split(clf, X_test, y_test )
    return clf, transformer


def classifyRealData(classifier, real_data, scalar):
    test_X = real_data.loc[:, ['height',                             
                              'mappingCountHeuristic',                             
                              'lostGene2Heuristic',
                              'lostGene1Heuristic',
                              'lostGene3Heuristic',                             
                              'geneFrequencyHeuristic'
                              ]]


    test_X[['height']] = scalar.fit_transform(test_X[['height']])

    output = classifier.predict(test_X)
    additive_count = 0
    additive_distance = 0
    replacing_count = 0
    replacing_distance = 0
    currentHeuristicAdditiveCount = 0
    currentHeuristicAdditiveDistance = 0
    currentHeuristicReplacingCount = 0
    currentHeuristicReplacingDistance = 0
    lossHeuristicAdditiveCount = 0
    lossHeuristicAdditiveDistance = 0
    lossHeuristicReplacingCount = 0
    lossHeuristicReplacingDistance = 0

    currentHeuristicLabel = real_data.loc[:, "currentHeuristicLabel"]
    lossHeuristicLabel = real_data.loc[:, "lossHeuristic"]

    for i in range(len(output)):
        if output[i]==1:
            additive_count +=1
            additive_distance+= real_data.iloc[i]['distanceOfTransfer']

        elif output[i]==0:
            replacing_count +=1
            replacing_distance += real_data.iloc[i]['distanceOfTransfer']

        if currentHeuristicLabel[i]==1:
            currentHeuristicAdditiveCount+=1
            currentHeuristicAdditiveDistance+= real_data.iloc[i]['distanceOfTransfer']

        elif currentHeuristicLabel[i]==0:
             currentHeuristicReplacingCount+=1
             currentHeuristicReplacingDistance+= real_data.iloc[i]['distanceOfTransfer']

        if lossHeuristicLabel[i]==1:
            lossHeuristicAdditiveCount+=1
            lossHeuristicAdditiveDistance+= real_data.iloc[i]['distanceOfTransfer']

        elif lossHeuristicLabel[i]==0:
            lossHeuristicReplacingCount+=1
            lossHeuristicReplacingDistance+= real_data.iloc[i]['distanceOfTransfer']


    print("Average additive by ML method: " + str(additive_count/len(output)*100.0))
    print("Average replacing by ML method: " + str(replacing_count / len(output) * 100.0))

    print("Average additive distance by ML method: " + str(additive_distance*1.0 / (additive_count)) )
    print("Average replacing distance by ML method: " + str(replacing_distance *1.0/ (replacing_count) ))

    print()
    print("Average additive by current heuristic: " + str(currentHeuristicAdditiveCount / len(output) * 100.0))
    print("Average replacing by current heuristic: " + str(currentHeuristicReplacingCount / len(output) * 100.0))
    print("Average additive distance by current heuristic: " + str(currentHeuristicAdditiveDistance * 1.0 / (currentHeuristicAdditiveCount)))
    print("Average replacing distance by current heuristic: " + str(currentHeuristicReplacingDistance * 1.0 / (currentHeuristicReplacingCount)))

    print()
    print("Average additive by loss heuristic: " + str(lossHeuristicAdditiveCount / len(output) * 100.0))
    print("Average replacing by loss heuristic: " + str(lossHeuristicReplacingCount / len(output) * 100.0))
    print("Average additive distance by loss heuristic: " + str(lossHeuristicAdditiveDistance * 1.0 / (lossHeuristicAdditiveCount)))
    print("Average replacing distance by loss heuristic: " + str(lossHeuristicReplacingDistance * 1.0 / (lossHeuristicReplacingCount)))
def classifySingleIO(classifier,test_data,  scalar):
    test_X = test_data.loc[:, ['height',                             
                              'mappingCountHeuristic',                             
                              'lostGene2Heuristic',
                              'lostGene1Heuristic',
                              'lostGene3Heuristic',                             
                              'geneFrequencyHeuristic'
                              ]]
    test_Y = test_data.trueLabel



    test_X[['height']] = scalar.fit_transform(test_X[['height']])

    output = classifier.predict(test_X)
    outputString = "";
    for val in output:
        outputString+=str(val)
    #print(len(output))
    print(outputString)

def removeFiles():
    if os.path.isfile("ML_output.txt"):
        os.remove("ML_output.txt")
    if os.path.isfile("current_heuristic_output.txt"):
        os.remove("current_heuristic_output.txt")
    if os.path.isfile("loss_heuristic_output.txt"):
        os.remove("loss_heuristic_output.txt")
    if os.path.isfile("ML_ranger_dependant_output.txt"):
        os.remove("ML_ranger_dependant_output.txt")
    if os.path.isfile("current_heuristic_ranger_dependant_output.txt"):
        os.remove("current_heuristic_ranger_dependant_output.txt")
    if os.path.isfile("loss_heuristic_ranger_dependant_output.txt"):
        os.remove("loss_heuristic_ranger_dependant_output.txt")
    

rand = random.randint(0, 9)
'''
low_replacing_test_data = pd.read_csv('low_replacing_test.txt', index_col='id', true_values='trueLabel')
low_mixed_test_data = pd.read_csv('low_mixed_test.txt', index_col='id', true_values='trueLabel')
low_additive_test_data = pd.read_csv('low_additive_test.txt', index_col='id', true_values='trueLabel')
medium_additive_test_data = pd.read_csv('medium_additive_test.txt', index_col='id', true_values='trueLabel')
medium_mixed_test_data = pd.read_csv('medium_mixed_test.txt', index_col='id', true_values='trueLabel')
medium_replacing_test_data = pd.read_csv('medium_replacing_test.txt', index_col='id', true_values='trueLabel')
high_additive_test_data = pd.read_csv('high_additive_test.txt', index_col='id', true_values='trueLabel')
high_mixed_test_data = pd.read_csv('high_mixed_test.txt', index_col='id', true_values='trueLabel')
real_data = pd.read_csv('realData2Across.txt', index_col='id', true_values='trueLabel')
high_replacing_test_data = pd.read_csv('high_replacing_test.txt', index_col='id', true_values='trueLabel')
'''
training_data = pd.read_csv('all_training'  + '.txt', index_col='id', true_values='trueLabel')
if path.exists("saved_classifier.joblib"):
    classifier = load('saved_classifier.joblib')
    scalar = preprocessing.MinMaxScaler()
else:
    classifier, scalar = get_classifier(training_data)
dump(classifier, 'saved_classifier.joblib')
#print(classifier.feature_importances_)


removeFiles()
''' 
all_test_data = [low_additive_test_data,
                 low_mixed_test_data,
                 low_replacing_test_data,
                 medium_additive_test_data,
                 medium_mixed_test_data,
                 medium_replacing_test_data,
                 high_additive_test_data,
                 high_mixed_test_data,
                 high_replacing_test_data
                 ]

i = 1

for test_data in all_test_data:
    classify(classifier, test_data, i, scalar)
    i = i + 1
'''

single_IO_test_data = pd.read_csv(sys.argv[1], index_col='id', true_values='trueLabel')
classifySingleIO(classifier,single_IO_test_data,  scalar )
#classifyRealData(classifier,real_data, scalar)


if os.path.isfile(sys.argv[1]):
        os.remove(sys.argv[1])