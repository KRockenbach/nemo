# Environment Setup
Install [Mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html) to manage environments.

Install GNU:parallel:
```
apt install parallel
```

Setup environments:
```
cd ../dev_envs
mamba env create -f prep.yaml
# for GPU support
CONDA_OVERRIDE_CUDA="11.8" mamba env create -f ../envs/nemo.yaml
# for CPU support
mamba env create -f ../envs/nemo_cpu.yaml
mamba env create -f modisco.yaml
```


# Preprocessing:
`bash run_setup.sh`

# Hyperparameter optimization
`bash tune_hyperparams.sh 0 optimize_architecture.py architecture_optimized 1100 200 False`  
`CUDA_VISIBLE_DEVICES=0 python -m nemo.hyperparam_tuning.validate_architecture.py architectire_validated`  
`bash tune_hyperparams.sh 0 finetune_architecture.py architecture_finetuned 1100 200 False`  
`bash tune_hyperparams.sh 0 optimize_learning.py learning_optimized 510 300 True`  
`bash tune_hyperparams.sh 0 finetune_architecture.py learning_finetuned 510 300 True`  
`bash tune_hyperparams.sh 0 final_touches.py nemo 310 310 True`  
`CUDA_VISIBLE_DEVICES=0 python -m nemo.hyperparam_tuning.validate_final_architecture final_architectire_validated`

# Model training
`bash train_Bnapus.sh`  
`bash train_Athaliana.sh`

# Make predictions and group genes
`bash eval.sh`

# Attributions
`bash attribution.sh`

# Group importance
`bash importance.sh`

# In silico mutation analysis
`bash motif_insertion_pipeline.sh`

# Validation on real world data
`bash validation.sh`
