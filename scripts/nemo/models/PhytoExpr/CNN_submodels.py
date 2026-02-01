# Original author: Hai Wang (https://doi.org/10.1073/pnas.2319811121), modified by Kevin Rockenbach
# Original files: https://figshare.com/articles/dataset/code_zip/24417076?file=42825613
# License: https://creativecommons.org/licenses/by/4.0/

from tensorflow.keras.layers import *
from tensorflow.keras.models import Model
from tensorflow.keras.activations import relu


def conv_block(x,filters,kernel_size,padding,dilation_rate,pool_size,dropout):
   # changed activation from relu to leaky relu (alpha = 0.1)
   x=Conv1D(int(float(filters)),kernel_size=int(float(kernel_size)),padding=padding,dilation_rate=int(float(dilation_rate)))(x)
   x=relu(x, alpha=0.1)
   x=Conv1D(int(float(filters)),kernel_size=int(float(kernel_size)),padding=padding,dilation_rate=int(float(dilation_rate)))(x)
   x=relu(x, alpha=0.1)
   x=MaxPooling1D(pool_size=int(float(pool_size)),padding=padding)(x)
   x=Dropout(float(dropout))(x)
   return x

def dense_block(x,units,dropout):
   # changed activation from relu to leaky relu (alpha = 0.1)
   x=Dense(int(float(units)))(x)
   x=relu(x, alpha=0.1)
   x=Dropout(float(dropout))(x)
   return x

def dense_block_end(x,units,dropout,activation,name):
   x=Dense(int(float(units)),activation=activation,name=name)(x)
   return x

#################################################################################################
def build_C2D3(x,outside,parameters,prefix):
   [filtersC1,kernel_sizeC1,paddingC1,dilation_rateC1,pool_sizeC1,dropoutC1, \
   filtersC2,kernel_sizeC2,paddingC2,dilation_rateC2,pool_sizeC2,dropoutC2, \
   _,_,_,_,_,_, \
   _,_,_,_,_,_, \
   unitsB1D1,dropoutB1D1,unitsB1D2,dropoutB1D2]=parameters

   if outside == 3000:
       x = Cropping1D(1000)(x) 
   elif outside == 2000:
       x = Cropping1D(2000)(x) 
   else:
       x = Cropping1D(0)(x)

   x=conv_block(x,filtersC1,kernel_sizeC1,paddingC1,dilation_rateC1,pool_sizeC1,dropoutC1)
   x=conv_block(x,filtersC2,kernel_sizeC2,paddingC2,dilation_rateC2,pool_sizeC2,dropoutC2)

   x=Flatten(name=f'{prefix}flatten')(x)

   x=dense_block(x,unitsB1D1,dropoutB1D1)
   x=dense_block(x,unitsB1D2,dropoutB1D2)
   # output activation changed from relu to linear
   x=dense_block_end(x,1,0,activation='linear',name=f'{prefix}tpm_output')
   return x


def build_C3D3(x,outside,parameters,prefix):
   [filtersC1,kernel_sizeC1,paddingC1,dilation_rateC1,pool_sizeC1,dropoutC1, \
   filtersC2,kernel_sizeC2,paddingC2,dilation_rateC2,pool_sizeC2,dropoutC2, \
   filtersC3,kernel_sizeC3,paddingC3,dilation_rateC3,pool_sizeC3,dropoutC3, \
   _,_,_,_,_,_, \
   unitsB1D1,dropoutB1D1,unitsB1D2,dropoutB1D2]=parameters

   if outside == 3000:
       x = Cropping1D(1000)(x) 
   elif outside == 2000:
       x = Cropping1D(2000)(x) 
   else:
       x = Cropping1D(0)(x)

   x=conv_block(x,filtersC1,kernel_sizeC1,paddingC1,dilation_rateC1,pool_sizeC1,dropoutC1)
   x=conv_block(x,filtersC2,kernel_sizeC2,paddingC2,dilation_rateC2,pool_sizeC2,dropoutC2)
   x=conv_block(x,filtersC3,kernel_sizeC3,paddingC3,dilation_rateC3,pool_sizeC3,dropoutC3)

   x=Flatten(name=f'{prefix}flatten')(x)

   x=dense_block(x,unitsB1D1,dropoutB1D1)
   x=dense_block(x,unitsB1D2,dropoutB1D2)
   # output activation changed from relu to linear
   x=dense_block_end(x,1,0,activation='linear',name=f'{prefix}tpm_output')
   return x


def build_C4D3(x,outside,parameters,prefix):
   [filtersC1,kernel_sizeC1,paddingC1,dilation_rateC1,pool_sizeC1,dropoutC1, \
   filtersC2,kernel_sizeC2,paddingC2,dilation_rateC2,pool_sizeC2,dropoutC2, \
   filtersC3,kernel_sizeC3,paddingC3,dilation_rateC3,pool_sizeC3,dropoutC3, \
   filtersC4,kernel_sizeC4,paddingC4,dilation_rateC4,pool_sizeC4,dropoutC4, \
   unitsB1D1,dropoutB1D1,unitsB1D2,dropoutB1D2]=parameters

   if outside == 3000:
       x = Cropping1D(1000)(x) 
   elif outside == 2000:
       x = Cropping1D(2000)(x) 
   else:
       x = Cropping1D(0)(x)

   x=conv_block(x,filtersC1,kernel_sizeC1,paddingC1,dilation_rateC1,pool_sizeC1,dropoutC1)
   x=conv_block(x,filtersC2,kernel_sizeC2,paddingC2,dilation_rateC2,pool_sizeC2,dropoutC2)
   x=conv_block(x,filtersC3,kernel_sizeC3,paddingC3,dilation_rateC3,pool_sizeC3,dropoutC3)
   x=conv_block(x,filtersC4,kernel_sizeC4,paddingC4,dilation_rateC4,pool_sizeC4,dropoutC4)

   x=Flatten(name=f'{prefix}flatten')(x)

   x=dense_block(x,unitsB1D1,dropoutB1D1)
   x=dense_block(x,unitsB1D2,dropoutB1D2)
   # output activation changed from relu to linear
   x=dense_block_end(x,1,0,activation='linear',name=f'{prefix}tpm_output')
   return x


#################################################################################################
def build_C2D2(x,outside,parameters,prefix):
   [filtersC1,kernel_sizeC1,paddingC1,dilation_rateC1,pool_sizeC1,dropoutC1, \
   filtersC2,kernel_sizeC2,paddingC2,dilation_rateC2,pool_sizeC2,dropoutC2, \
   _,_,_,_,_,_, \
   _,_,_,_,_,_, \
   unitsB1D1,dropoutB1D1,_,_]=parameters

   if outside == 3000:
       x = Cropping1D(1000)(x) 
   elif outside == 2000:
       x = Cropping1D(2000)(x) 
   else:
       x = Cropping1D(0)(x)

   x=conv_block(x,filtersC1,kernel_sizeC1,paddingC1,dilation_rateC1,pool_sizeC1,dropoutC1)
   x=conv_block(x,filtersC2,kernel_sizeC2,paddingC2,dilation_rateC2,pool_sizeC2,dropoutC2)

   x=Flatten(name=f'{prefix}flatten')(x) 

   x=dense_block(x,unitsB1D1,dropoutB1D1)
   # output activation changed from relu to linear
   x=dense_block_end(x,1,0,activation='linear',name=f'{prefix}tpm_output')
   return x


def build_C3D2(x,outside,parameters,prefix):
   [filtersC1,kernel_sizeC1,paddingC1,dilation_rateC1,pool_sizeC1,dropoutC1, \
   filtersC2,kernel_sizeC2,paddingC2,dilation_rateC2,pool_sizeC2,dropoutC2, \
   filtersC3,kernel_sizeC3,paddingC3,dilation_rateC3,pool_sizeC3,dropoutC3, \
   _,_,_,_,_,_, \
   unitsB1D1,dropoutB1D1,_,_]=parameters

   if outside == 3000:
       x = Cropping1D(1000)(x) 
   elif outside == 2000:
       x = Cropping1D(2000)(x) 
   else:
       x = Cropping1D(0)(x)

   x=conv_block(x,filtersC1,kernel_sizeC1,paddingC1,dilation_rateC1,pool_sizeC1,dropoutC1)
   x=conv_block(x,filtersC2,kernel_sizeC2,paddingC2,dilation_rateC2,pool_sizeC2,dropoutC2)
   x=conv_block(x,filtersC3,kernel_sizeC3,paddingC3,dilation_rateC3,pool_sizeC3,dropoutC3)

   x=Flatten(name=f'{prefix}flatten')(x)

   x=dense_block(x,unitsB1D1,dropoutB1D1)
   # output activation changed from relu to linear
   x=dense_block_end(x,1,0,activation='linear',name=f'{prefix}tpm_output')
   return x


def build_C4D2(x,outside,parameters,prefix):
   [filtersC1,kernel_sizeC1,paddingC1,dilation_rateC1,pool_sizeC1,dropoutC1, \
   filtersC2,kernel_sizeC2,paddingC2,dilation_rateC2,pool_sizeC2,dropoutC2, \
   filtersC3,kernel_sizeC3,paddingC3,dilation_rateC3,pool_sizeC3,dropoutC3, \
   filtersC4,kernel_sizeC4,paddingC4,dilation_rateC4,pool_sizeC4,dropoutC4, \
   unitsB1D1,dropoutB1D1,_,_]=parameters

   if outside == 3000:
       x = Cropping1D(1000)(x) 
   elif outside == 2000:
       x = Cropping1D(2000)(x) 
   else:
       x = Cropping1D(0)(x)

   x=conv_block(x,filtersC1,kernel_sizeC1,paddingC1,dilation_rateC1,pool_sizeC1,dropoutC1)  
   x=conv_block(x,filtersC2,kernel_sizeC2,paddingC2,dilation_rateC2,pool_sizeC2,dropoutC2)  
   x=conv_block(x,filtersC3,kernel_sizeC3,paddingC3,dilation_rateC3,pool_sizeC3,dropoutC3)  
   x=conv_block(x,filtersC4,kernel_sizeC4,paddingC4,dilation_rateC4,pool_sizeC4,dropoutC4)  

   x=Flatten(name=f'{prefix}flatten')(x)

   x=dense_block(x,unitsB1D1,dropoutB1D1)
   # output activation changed from relu to linear
   x=dense_block_end(x,1,0,activation='linear',name=f'{prefix}tpm_output')
   return x



