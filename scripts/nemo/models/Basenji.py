'''
MIT License

Copyright (c) 2025 liulifen

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
'''

import tensorflow as tf
import numpy as np
from tensorflow.keras.layers import Input, \
                                    BatchNormalization, \
                                    Conv1D, \
                                    Dropout, \
                                    MaxPooling1D
from tensorflow.keras.models import Model,load_model

# best found input configuration for Arabidopsis: 2.5K upstream of TSS and 2.5K downstream of TSS!

channels_num = 512

def exponential_linspace_int(initial_value, target_value, num_layers):
    factor = (target_value / initial_value) ** (1 / num_layers)
    values = []
    # Calculate and store values
    value = initial_value
    for _ in range(num_layers+1):
        values.append(np.round(value))
        value *= factor
    return values[1:]


class GELU(tf.keras.layers.Layer):
    def __init__(self, **kwargs):
        super(GELU, self).__init__(**kwargs)
    def __call__(self, x):
       # return tf.keras.activations.sigmoid(1.702 * x) * x
        return tf.keras.activations.sigmoid(tf.constant(1.702) * x) * x

def conv_block(x, C=channels_num, W=1, D=1, kernel_initializer='he_normal', l2_scale=0):
    x = Conv1D(
        filters=int(C),
        kernel_size=W,
        #strides=strides,
        padding='same',
        dilation_rate=D,
        kernel_initializer=kernel_initializer,
        kernel_regularizer=tf.keras.regularizers.l2(l2_scale))(x)
    x = BatchNormalization()(x)
    x = GELU()(x)

    return x

def residual_block(input, D=1):
    x = conv_block(input, C=0.5*channels_num, W=3, D=D)
    x = conv_block(x, C=channels_num, W=1, D=1)
    x = Dropout(0.3)(x)
    output = x + input
    return output

def basenji_model(input_shape=[5000, 4], W = 15,L=11):
    inputs = Input(input_shape)
    x = Conv1D(
        filters=int(0.375*channels_num),
        kernel_size=W,
        padding='same',
        activation='exponential',
        dilation_rate=1,
        kernel_initializer='he_normal',
        kernel_regularizer=tf.keras.regularizers.l2(0))(inputs)
    x = BatchNormalization()(x)
    x = GELU()(x)
    x = MaxPooling1D(pool_size= 3 )(x)  
    Ci_steps = exponential_linspace_int(0.5*channels_num, channels_num, 6) 
    for Ci in Ci_steps:
        x = conv_block(x=x, C=Ci, W=5, D=1)
        x = MaxPooling1D(pool_size=2)(x)
    
    Di=[1,2,3,4]
    for i in range(len(Di)): 
        x = residual_block(x, D=Di[i]) 


    x = conv_block(x, C=channels_num//2, W=1, D=1)
    x = Dropout(0.05)(x)
    x = GELU()(x)
    x = Conv1D(filters = 1,
               kernel_size=1, 
               padding='same')(x)
    x = tf.keras.layers.Flatten()(x)
    x= tf.keras.layers.Dense(1)(x)
    model = Model(inputs = inputs, outputs = x)
    return model
