"""
This script is used for ResNet model trainning.

"""

from cgi import test
from torchvision import models
import torch.nn as nn
import torch.optim as optim
from IPython.display import clear_output
import math
import torch
import os


def training(device,trainloader,testloader,evaluate,train,learning_rate,num_epochs, save_checkpoint, model_name):
    
    save_path = os.path.join('Models', model_name + ".pt")

    model= models.resnet50(pretrained=True)

    num_ftrs = model.fc.in_features
    model.fc = nn.Linear(num_ftrs, 2)
    model = model.to(device)

    start_epoch = 0

    optimizer = optim.Adam(model.parameters(), lr = learning_rate)
    loss_fun = nn.CrossEntropyLoss()

    training_loss_logger = []
    # validation_loss_logger = []

    training_acc_logger = []
    # validation_acc_logger = []
    testing_acc_logger = []
    testing_loss_logger = []

    highest_test_acc = -math.inf
    valid_acc = 0   # Just setting a random acc to store in model, since we do not have validation images yet

    #This cell implements our training loop
    for epoch in range(start_epoch, num_epochs):
        
        #call the training function and pass training dataloader etc
        training_loss_logger = train(model, device, trainloader, optimizer, loss_fun, training_loss_logger)
        
        #call the evaluate function and pass the dataloader for both ailidation and training
        train_acc, _ = evaluate(model, device, trainloader, loss_fun)
        training_acc_logger.append(train_acc)
        
    #     valid_acc, valid_loss = evaluate(model, device, validloader, loss_fun)
    #     validation_loss_logger.append(valid_loss)
    #     validation_acc_logger.append(valid_acc)
        
        test_acc, test_loss = evaluate(model, device, testloader, loss_fun)
        testing_loss_logger.append(test_loss)
        testing_acc_logger.append(test_acc)

        #If this model has the highest performace on the validation set 
        #then save a checkpoint
        #{} define a dictionary, each entry of the dictionary is indexed with a string
        if (test_acc >= highest_test_acc):
            highest_test_acc = test_acc
            highest_test_acc_epoch = epoch

            if save_checkpoint:
                print("Saving Model")
                torch.save({
                    'epoch':                 epoch,
                    'model_state_dict':      model.state_dict(),
                    'optimizer_state_dict':  optimizer.state_dict(), 
                    'train_acc':             train_acc,
                    'valid_acc':             valid_acc,
                }, save_path)
        
        clear_output(True)
        print(f'| Epoch: {epoch+1:02} | Train Acc: {train_acc:05.2f}% | Val. Acc: {valid_acc:05.2f}% | Test Acc: {test_acc:05.2f}%')

    print("Training Complete")
    print('The final test accuracy is ',test_acc)
    print('The highest test accuracy is ',highest_test_acc)


    return training_loss_logger, training_acc_logger, testing_acc_logger, testing_loss_logger, highest_test_acc,model