# Original author: Hai Wang (https://doi.org/10.1073/pnas.2319811121), modified by Kevin Rockenbach
# Original files: https://figshare.com/articles/dataset/code_zip/24417076?file=42825613
# License: https://creativecommons.org/licenses/by/4.0/

from keras.layers import *
from keras.models import Model

def conv_block(x,filters,kernel_size,padding,dilation_rate,pool_size,dropout):
   x=Conv1D(int(float(filters)),kernel_size=int(float(kernel_size)),padding=padding,dilation_rate=int(float(dilation_rate)),activation='relu')(x)
   x=Conv1D(int(float(filters)),kernel_size=int(float(kernel_size)),padding=padding,dilation_rate=int(float(dilation_rate)),activation='relu')(x)
   x=MaxPooling1D(pool_size=int(float(pool_size)),padding=padding)(x)
   x=Dropout(float(dropout))(x)
   return x

def dense_block(x,units,dropout):
   x=Dense(int(float(units)),activation='relu')(x)
   x=Dropout(float(dropout))(x)
   return x

def dense_block_end(x,units,dropout,activation,name):
   x=Dense(int(float(units)),activation=activation,name=name)(x)
   return x

#################################################################################################
def build_C2D3(x,parameters,prefix):
   [filtersC1,kernel_sizeC1,paddingC1,dilation_rateC1,pool_sizeC1,dropoutC1, \
   filtersC2,kernel_sizeC2,paddingC2,dilation_rateC2,pool_sizeC2,dropoutC2, \
   _,_,_,_,_,_, \
   _,_,_,_,_,_, \
   unitsB1D1,dropoutB1D1,unitsB1D2,dropoutB1D2]=parameters

   x=conv_block(x,filtersC1,kernel_sizeC1,paddingC1,dilation_rateC1,pool_sizeC1,dropoutC1)
   x=conv_block(x,filtersC2,kernel_sizeC2,paddingC2,dilation_rateC2,pool_sizeC2,dropoutC2)

   x=Flatten(name='flatten')(x)

   x=dense_block(x,unitsB1D1,dropoutB1D1)
   x=dense_block(x,unitsB1D2,dropoutB1D2)
   x=dense_block_end(x,1,0,activation='relu',name=f'{prefix}tpm_output')
   return x


def build_C3D3(x,parameters,prefix):
   [filtersC1,kernel_sizeC1,paddingC1,dilation_rateC1,pool_sizeC1,dropoutC1, \
   filtersC2,kernel_sizeC2,paddingC2,dilation_rateC2,pool_sizeC2,dropoutC2, \
   filtersC3,kernel_sizeC3,paddingC3,dilation_rateC3,pool_sizeC3,dropoutC3, \
   _,_,_,_,_,_, \
   unitsB1D1,dropoutB1D1,unitsB1D2,dropoutB1D2]=parameters

   x=conv_block(x,filtersC1,kernel_sizeC1,paddingC1,dilation_rateC1,pool_sizeC1,dropoutC1)
   x=conv_block(x,filtersC2,kernel_sizeC2,paddingC2,dilation_rateC2,pool_sizeC2,dropoutC2)
   x=conv_block(x,filtersC3,kernel_sizeC3,paddingC3,dilation_rateC3,pool_sizeC3,dropoutC3)

   x=Flatten(name='flatten')(x)

   x=dense_block(x,unitsB1D1,dropoutB1D1)
   x=dense_block(x,unitsB1D2,dropoutB1D2)
   x=dense_block_end(x,1,0,activation='relu',name=f'{prefix}tpm_output')
   return x


def build_C4D3(x,parameters,prefix):
   [filtersC1,kernel_sizeC1,paddingC1,dilation_rateC1,pool_sizeC1,dropoutC1, \
   filtersC2,kernel_sizeC2,paddingC2,dilation_rateC2,pool_sizeC2,dropoutC2, \
   filtersC3,kernel_sizeC3,paddingC3,dilation_rateC3,pool_sizeC3,dropoutC3, \
   filtersC4,kernel_sizeC4,paddingC4,dilation_rateC4,pool_sizeC4,dropoutC4, \
   unitsB1D1,dropoutB1D1,unitsB1D2,dropoutB1D2]=parameters

   x=conv_block(x,filtersC1,kernel_sizeC1,paddingC1,dilation_rateC1,pool_sizeC1,dropoutC1)
   x=conv_block(x,filtersC2,kernel_sizeC2,paddingC2,dilation_rateC2,pool_sizeC2,dropoutC2)
   x=conv_block(x,filtersC3,kernel_sizeC3,paddingC3,dilation_rateC3,pool_sizeC3,dropoutC3)
   x=conv_block(x,filtersC4,kernel_sizeC4,paddingC4,dilation_rateC4,pool_sizeC4,dropoutC4)

   x=Flatten(name='flatten')(x)

   x=dense_block(x,unitsB1D1,dropoutB1D1)
   x=dense_block(x,unitsB1D2,dropoutB1D2)
   x=dense_block_end(x,1,0,activation='relu',name=f'{prefix}tpm_output')
   return x


#################################################################################################
def build_C2D2(x,parameters,prefix):
   [filtersC1,kernel_sizeC1,paddingC1,dilation_rateC1,pool_sizeC1,dropoutC1, \
   filtersC2,kernel_sizeC2,paddingC2,dilation_rateC2,pool_sizeC2,dropoutC2, \
   _,_,_,_,_,_, \
   _,_,_,_,_,_, \
   unitsB1D1,dropoutB1D1,_,_]=parameters

   x=conv_block(x,filtersC1,kernel_sizeC1,paddingC1,dilation_rateC1,pool_sizeC1,dropoutC1)
   x=conv_block(x,filtersC2,kernel_sizeC2,paddingC2,dilation_rateC2,pool_sizeC2,dropoutC2)

   x=Flatten(name='flatten')(x) 

   x=dense_block(x,unitsB1D1,dropoutB1D1)
   x=dense_block_end(x,1,0,activation='relu',name=f'{prefix}tpm_output')
   return x


def build_C3D2(x,parameters,prefix):
   [filtersC1,kernel_sizeC1,paddingC1,dilation_rateC1,pool_sizeC1,dropoutC1, \
   filtersC2,kernel_sizeC2,paddingC2,dilation_rateC2,pool_sizeC2,dropoutC2, \
   filtersC3,kernel_sizeC3,paddingC3,dilation_rateC3,pool_sizeC3,dropoutC3, \
   _,_,_,_,_,_, \
   unitsB1D1,dropoutB1D1,_,_]=parameters

   x=conv_block(x,filtersC1,kernel_sizeC1,paddingC1,dilation_rateC1,pool_sizeC1,dropoutC1)
   x=conv_block(x,filtersC2,kernel_sizeC2,paddingC2,dilation_rateC2,pool_sizeC2,dropoutC2)
   x=conv_block(x,filtersC3,kernel_sizeC3,paddingC3,dilation_rateC3,pool_sizeC3,dropoutC3)

   x=Flatten(name='flatten')(x)

   x=dense_block(x,unitsB1D1,dropoutB1D1)
   x=dense_block_end(x,1,0,activation='relu',name=f'{prefix}tpm_output')
   return x


def build_C4D2(x,parameters,prefix):
   [filtersC1,kernel_sizeC1,paddingC1,dilation_rateC1,pool_sizeC1,dropoutC1, \
   filtersC2,kernel_sizeC2,paddingC2,dilation_rateC2,pool_sizeC2,dropoutC2, \
   filtersC3,kernel_sizeC3,paddingC3,dilation_rateC3,pool_sizeC3,dropoutC3, \
   filtersC4,kernel_sizeC4,paddingC4,dilation_rateC4,pool_sizeC4,dropoutC4, \
   unitsB1D1,dropoutB1D1,_,_]=parameters

   x=conv_block(x,filtersC1,kernel_sizeC1,paddingC1,dilation_rateC1,pool_sizeC1,dropoutC1)  
   x=conv_block(x,filtersC2,kernel_sizeC2,paddingC2,dilation_rateC2,pool_sizeC2,dropoutC2)  
   x=conv_block(x,filtersC3,kernel_sizeC3,paddingC3,dilation_rateC3,pool_sizeC3,dropoutC3)  
   x=conv_block(x,filtersC4,kernel_sizeC4,paddingC4,dilation_rateC4,pool_sizeC4,dropoutC4)  

   x=Flatten(name='flatten')(x)

   x=dense_block(x,unitsB1D1,dropoutB1D1)
   x=dense_block_end(x,1,0,activation='relu',name=f'{prefix}tpm_output')
   return x



