from numpy import random, array, zeros, empty
import os
from sklearn.model_selection import KFold
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from xgboost import XGBClassifier
import joblib


def readTest(filename):
	data = open(filename, "r")
	lines = data.readlines()

	data_feature = []

	for line in lines:
		line = line.strip()
		elements = line.split("\t")
		features = [float(elements[i]) for i in range(0, len(elements))]
		data_feature.append(features)

	data_feature = array(data_feature)

	return data_feature


def getRange(feature):
	feature_range = zeros((feature.shape[1], 2))

	for i in range(0, feature.shape[1]):
		col = feature[:, i]
		min_val = float(min(col))
		max_val = float(max(col))
		feature_range[i, :] = [min_val, max_val]

	return feature_range


def svm_scale(feature, ranges, lower, upper):
	feature_scaled = zeros((feature.shape[0], feature.shape[1]))

	for i in range(0, feature.shape[1]-1):
		feat = feature[:, i]
		min_val = ranges[i, 0]
		max_val = ranges[i, 1]
		if ((max_val - min_val) != 0):
			feat_std = (feat - min_val) / (max_val - min_val)
			feature_scaled[:, i] = feat_std * (upper - lower) + lower
		else:
			feature_scaled[:, i] = zeros((feature.shape[0],))

	return feature_scaled


def stacking_predict(features_test):
	ntest = features_test.shape[0]
	second_test = zeros((ntest, 2))
	svm_test_kf = empty((5, ntest))
	boost_test_kf = empty((5, ntest))

	train_range = joblib.load(cwd + '/models/train_range.joblib')

	features_test_scale = svm_scale(features_test, train_range, 0, 1)

	for i in range(0,5):
		filename1 = cwd + '/models/SVM_model{}.joblib'.format(i)
		svm_classifier = joblib.load(filename1)
		svm_test_kf[i,:] = svm_classifier.predict_proba(features_test_scale)[:, 1]

	svm_test = svm_test_kf.mean(axis = 0)
	second_test[:, 0] = svm_test

	for j in range(0,5):
		filename2 = cwd + '/models/XGBoost_model{}.joblib'.format(j)
		boost_classifier = joblib.load(filename2)
		boost_test_kf[j, :] = boost_classifier.predict_proba(features_test)[:, 1]

	boost_test = boost_test_kf.mean(axis = 0)
	second_test[:, 1] = boost_test

	second_filename = cwd + '/models/second_logistic_model.joblib'
	second_layer_classifier = joblib.load(second_filename)
	label_pred = second_layer_classifier.predict(second_test)
	prob_pred = second_layer_classifier.predict_proba(second_test)[:, 1] # row logistic regression prob [0.0, 1.0]

	return label_pred, prob_pred

temp_dir = os.environ["SGDESIGNER_TEMP_DIR"]
cwd = os.getcwd()

# customer's dataset
test_file = temp_dir + "/custom_features_v2.0.txt"
# result file
result_file = temp_dir + "/custom_prediction_result_v2.0.txt"

test_feature = readTest(test_file)

label_pred, prob_pred = stacking_predict(test_feature)
prob_pred = [int(i*100+0.5) for i in prob_pred] # scale to 0-100 int with rounding
prob_pred = [((j - 30.5)/(76 - 30.5) * 100) for j in prob_pred] # re-normalize against empirical training distribution
prob_pred = [int(k +0.5) for k in prob_pred] # round to final integer


f = open(result_file, 'w')
try:
	f.write("labels\tprobabilities\n")
	for i, (label) in enumerate(label_pred):
		f.write("{}\t{}\n".format(label, prob_pred[i]))

finally:
	f.close()

