import sys, h5py, os
import pandas as pd
import tensorflow as tf
from tensorflow.keras.models import load_model, Model
from tensorflow.keras.layers import Concatenate, Dense, Input
from ..utils.model_utils import *


for modelname in ["Basenji-5K", "xpresso", "xpresso_no_halflife", "nemo", "PhytoExpr_CNN", "PhytoExpr_transformer"]:
    modeldir = os.path.join("..", "model_configs", modelname)

    if modelname == "Basenji-5K":
        from ..models.Basenji import basenji_model
        model = basenji_model()
    elif modelname == "xpresso" or modelname == "xpresso_no_halflife":
        from ..models.Xpresso import build_xpresso, build_xpresso_no_halflife
        if modelname == "xpresso":
            model = build_xpresso()
        else:
            model = build_xpresso_no_halflife()
    elif modelname == "nemo":
        from ..models.nemo import build_nemo
        model = build_nemo()
    elif modelname == "PhytoExpr_CNN":
        from ..models.PhytoExpr import CNN_submodels as sm
        cfg = pd.read_csv(os.path.join(modeldir, "ensemble_model_cfg.csv"), delimiter=";", header=0)
        model_types = cfg.loc[:,"type"].tolist()
        outside_intervals = cfg.loc[:,"outside_interval"].tolist()
        inside_intervals = cfg.loc[:,"inside_interval"].tolist()
        param_df = cfg.iloc[:,4:-1]
        x = []
        submodel_dict={}
        input_length=10000
        input=Input(shape=(input_length,4))
        for submodel_index in range(27):
            prefix = f"submodel{submodel_index}_"
            model_type = model_types[submodel_index] # C2D3, C3D3, C4D3, C2D2, C3D2, C4D2
            outside = outside_intervals[submodel_index]
            params = param_df.iloc[submodel_index,:].tolist()
            # build submodel
            exec(f"x{submodel_index} = sm.build_{model_type}(input,outside,params,prefix)")
            exec(f"x.append(x{submodel_index})")
        x=Concatenate()(x)
        x=Dense(1,name='second_layer_model_tpm')(x)
        model = Model(inputs=input,outputs=x)
    elif modelname == "PhytoExpr_transformer":
        from ..models.PhytoExpr.transformer import procheck
        input_length=10000
        input=Input(shape=(input_length, 4))
        Procheck=procheck(input)
        model=Model(inputs=input, outputs=Procheck)

    print(modelname + "\t" + str(model.count_params()))
