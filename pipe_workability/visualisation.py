""" 
This script is used to 
1. plot the loss and accuracy graphs for training, validation and test data
2. plot confusion matrix table 

"""

from matplotlib import pyplot as plt
import numpy as np
import pandas as pd

def plotting(epoch,training_loss_logger,training_acc_logger,testing_acc_logger,testing_loss_logger):
    """
    Function to plot loss and accuracy for traning and test data
    :param epoch: the number times that the learning algorithm will work through the entire training dataset
    :param *logger: a list contains loss or accuracy
    """
    # Training Loss
    plt.figure(figsize=(8, 4), dpi=80)
    plt.subplots_adjust(right = 1.3)
    plt.tight_layout()
    plt.subplot(1,2,1)
    train_x = np.linspace(0, epoch, len(training_loss_logger))
    plt.plot(train_x, training_loss_logger)
    plt.title('Training Loss Per Epoch')
    plt.xlabel('Epoch Number')
    plt.ylabel('Loss')

    # # Validation Loss
    # plt.subplot(1,2,2)
    # valid_x = np.linspace(0, epoch, len(validation_loss_logger))
    # plt.plot(valid_x, validation_loss_logger)
    # plt.title('Validation Loss Per Epoch')
    # plt.xlabel('Epoch Number')
    # plt.ylabel('Loss')

    # Testing Loss
    plt.subplot(1,2,2)
    test_x = np.linspace(0, epoch, len(testing_acc_logger))
    valid_x = np.linspace(0, epoch, len(testing_loss_logger))
    plt.plot(test_x, testing_loss_logger)
    plt.title('Validation Loss Per Epoch')
    plt.xlabel('Epoch Number')
    plt.ylabel('Loss')

    # Training Accuracy
    plt.figure(figsize=(11, 4), dpi=80)
    plt.subplots_adjust(right = 1.5)
    plt.tight_layout()
    plt.subplot(1,3,1)
    train_x = np.linspace(0, epoch, len(training_acc_logger))
    plt.plot(train_x, training_acc_logger)
    plt.title('Training Accuracy Per Epoch')
    plt.xlabel('Epoch Number')
    plt.ylabel('Accuracy')

    # # Validation Accuracy
    # plt.subplot(1,3,2)
    # valid_x = np.linspace(0, epoch, len(validation_acc_logger))
    # plt.plot(valid_x, validation_acc_logger)
    # plt.title('Validation Accuracy Per Epoch')
    # plt.xlabel('Epoch Number')
    # plt.ylabel('Accuracy')

    # Testing Accuracy
    plt.subplot(1,3,3)
    test_x = np.linspace(0, epoch, len(testing_acc_logger))
    plt.plot(test_x, testing_acc_logger)
    plt.title('Testing Accuracy Per Epoch')
    plt.xlabel('Epoch Number')
    plt.ylabel('Accuracy')


def confusion_matrix(cm):
    """
    Function to plot confusion matrix in table format
    :param cm: confusion matrix in array format
    :returns: confusion matrix and its corresponding sum in table format
    """
    # confusion matrix value calculation
    true_neg, false_pos, false_neg, true_pos = cm.ravel()

    # assign column and row labels
    col_labels = ('Predicted Positive', 'Predicted Negative', 'Sum')
    row_labels = ('Actual Positive', 'Actual Negative', 'Sum')

    array = [[true_pos,false_neg,true_pos+false_neg],
            [false_pos,true_neg,true_neg+false_pos],
            [true_pos+false_pos, true_neg+false_neg, true_pos+false_pos+true_neg+false_neg]]

    # create datafram
    df_cm = pd.DataFrame(array, index = [i for i in row_labels],
                            columns = [i for i in col_labels])
    
    return df_cm