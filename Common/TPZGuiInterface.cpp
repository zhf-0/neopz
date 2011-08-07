/**
 * @file
 * @brief Contains the TPZGuiInterface methods.
 */
//$Id: TPZGuiInterface.cpp,v 1.3 2010-07-23 04:54:25 phil Exp $

#include "TPZGuiInterface.h"
#include "pzerror.h"

TPZGuiInterface::TPZGuiInterface(){
	this->fCanceled = false;
	this->fMessage = "";
	this->fProgressBarPos = 0;
	this->fProgressBarMaxPos = 0;
	this->fProgressBarMinPos = 0;
}

TPZGuiInterface::~TPZGuiInterface(){
	///nothing to be done
}

void TPZGuiInterface::UpdateCaption(){
	std::cout << fMessage.c_str() << "\n"
	<< "Progress bar = " << fProgressBarPos << "/" << fProgressBarMaxPos
	<< "\n";
}

void TPZGuiInterface::Start(){
	std::cout << "Starting execution\n";
}

void TPZGuiInterface::End(){
	std::cout << "Execution finished\n";
}

void TPZGuiInterface::ShowErrorMessage(std::string message){
	PZError << message.c_str() << "\n";
	DebugStop();
}


void TPZGuiInterface::SetKilled(){
	this->fCanceled = true;
}


bool TPZGuiInterface::AmIKilled(){
	return this->fCanceled;
}



