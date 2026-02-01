# encoding: utf-8

# Original author: Hai Wang (https://doi.org/10.1073/pnas.2319811121), modified by Kevin Rockenbach
# Original files: https://figshare.com/articles/dataset/code_zip/24417076?file=42825613
# License: https://creativecommons.org/licenses/by/4.0/

from tensorflow.keras.layers import *
import math

def exponential_linspace_int(start, end, num, divisible_by = 1):
   def _round(x):
      return int(round(x / divisible_by) * divisible_by)
   base = math.exp(math.log(end / start) / (num - 1))
   return [_round(start * base**i) for i in range(num)]

def pooling_module(pool_size=6, kind='max'):
  if kind == 'max':
    return MaxPool1D(pool_size=pool_size, padding='same')
  else:
    raise ValueError('Invalid pooling type')

def BatchActivate(x):
    x = BatchNormalization(momentum=0.9, scale=True)(x)#momentum？
    x = Activation('gelu')(x)
    return x

def stem(x,Kernel_Size=15,Pad='same',Stride=1,Filter=128,pooling_type='max'):
    name='stem'
    x=Conv1D(filters=Filter,kernel_size=Kernel_Size,padding=Pad,strides=Stride, name=name)(x)
    x_forres=x
    x=convolution_block(x,Filter, size = 1, activation=True)
    x=residual_block(x,x_forres)
    x=pooling_module(pool_size=4, kind=pooling_type)(x)
    return x

# convblock
def convolution_block(x, filters, size=1, strides=1, padding='same', activation=True):
    if activation == True:
        x = BatchActivate(x)
    x = Conv1D(filters, size, strides=strides, padding=padding)(x)
    return x

# residual_block
def residual_block(x,short_x,batch_activate = False,conv_activate=True):
    x = Add()([short_x, x])
    if batch_activate:
        x = BatchActivate(x)
    return x

def conv_tower(x,tower_poolsize,conv_tower_num,dim,pooling_type='max'):
    filter=exponential_linspace_int(dim//2, dim, conv_tower_num, divisible_by = 128)
    for filter_num in filter:
        x_forres=convolution_block(x,filter_num, size = 5, activation=True,padding='same')
        x=convolution_block(x_forres,filter_num, size = 1, activation=True)
        x=residual_block(x,x_forres)
        x=pooling_module(pool_size=tower_poolsize, kind=pooling_type)(x)
        print(x.shape)
    return x

def mlp(x,dim,dropout_rate=0.3):
    x=LayerNormalization()(x)
    x=Dense(dim*2,activation='linear')(x)
    x=Dropout(0.15)(x)
    x=Activation('relu')(x)
    x=Dense(dim,activation='linear')(x)
    x=Dropout(dropout_rate)(x)
    return x

def MHA(x,dim,num_transformer_layers,heads=2,dropout_rate=0.3):
    for i in range(num_transformer_layers):
        x_forres=x
        x=LayerNormalization()(x)
        x=MultiHeadAttention(heads,32)(x,x)
        x=Dropout(dropout_rate)(x)
        x=residual_block(x,x_forres)
        x_forres=x
        x=mlp(x,dim,dropout_rate=dropout_rate)
        x=residual_block(x,x_forres)
    return x

def final(x,dim,crop_len=0,dropout_rate=0.3):
   x=Cropping1D(cropping=crop_len)(x)
   x=convolution_block(x,2*dim,size=1)
   x=Dropout(dropout_rate)(x)
   x=Activation('gelu')(x)
   x=Flatten()(x)
   x=Dense(1,activation = 'softplus',name='tpm_output')(x)
   return x

def procheck(x):
   stem_size = 15
   tower_poolsize = 4
   dim = 512
   conv_tower_num = 5
   num_transformer_layers = 6
   attention_heads = 10
   dropout_rate = 0.2

   dropout_rate=float(dropout_rate)
   Stem=stem(x,Kernel_Size=stem_size,Pad='same',Stride=1,Filter=dim//2)
   Conv_tower=conv_tower(Stem,tower_poolsize,conv_tower_num,dim)
   Mha=MHA(Conv_tower,dim,num_transformer_layers,heads=attention_heads,dropout_rate=dropout_rate)
   Final=final(Mha,dim,crop_len=0,dropout_rate=dropout_rate)
   return Final

