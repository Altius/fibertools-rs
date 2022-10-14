import xgboost as xgb
import numpy as np

train_csv = "models/temp/train.csv"
val_csv = "models/temp/val.csv"

read_data = False

if read_data:
    data_path = "models/temp/m6A_train_other_half_hifi.npz"
    train_val_data = np.load(data_path, allow_pickle=True)

    # Get the dictionary from the containing relevant data
    train_val_data = train_val_data["save_data_dict"][()]

    # Load training and validation features and labels
    X_train = train_val_data["X_train"]
    y_train = train_val_data["y_train"]
    X_val = train_val_data["X_val"]
    y_val = train_val_data["y_val"]

    # Reshape features to one dimension
    X_train = X_train.reshape(X_train.shape[0], X_train.shape[1] * X_train.shape[2])
    X_val = X_val.reshape(X_val.shape[0], X_val.shape[1] * X_val.shape[2])

    # Add labels before feature list
    y_train = y_train[:, 0]
    y_val = y_val[:, 0]
    final_train_data = np.concatenate([y_train[:, np.newaxis], X_train], axis=1)
    final_val_data = np.concatenate([y_val[:, np.newaxis], X_val], axis=1)
    print("saving")
    # Save data as CSV
    np.savetxt(train_csv, final_train_data, delimiter=",")
    np.savetxt(val_csv, final_val_data, delimiter=",")
    print("saved")

# Read in training and validation data
dtrain = xgb.DMatrix(f"{train_csv}?format=csv&label_column=0")
dval = xgb.DMatrix(f"{val_csv}?format=csv&label_column=0")

# Specify parameters via map
param = {"max_depth": 6, "objective": "binary:logistic", "eval_metric": "auc"}
num_round = 100
eval_s = [(dtrain, "train"), (dval, "validation")]

# Train model
bst = xgb.train(param, dtrain, num_round, evals=eval_s)

# Save model
bst.save_model("models/xgboost.0.81.bin")