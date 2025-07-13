import logging
from sklearn.metrics import classification_report, accuracy_score, confusion_matrix, mean_squared_error, mean_absolute_error
import pandas as pd
import datetime

# Configure logging: no timestamp in every line
logging.basicConfig(
    filename='results.log',
    level=logging.INFO,
    format='%(levelname)s - %(message)s'
)

def evaluate_and_log(model, X_test, y_test, sample_df):
    y_pred = model.predict(X_test)
    acc = accuracy_score(y_test, y_pred)
    report = classification_report(y_test, y_pred)
    cm = confusion_matrix(y_test, y_pred)

    mse = mean_squared_error(y_test, y_pred)
    mae = mean_absolute_error(y_test, y_pred)

    logging.info("--------- Evaluation Results ---------")
    logging.info("Timestamp: %s", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

    logging.info(f"Accuracy: {acc:.4f}")
    logging.info(f"Mean Squared Error (MSE): {mse:.4f}")
    logging.info(f"Mean Absolute Error (MAE): {mae:.4f}")
    logging.info("Classification Report:\n%s", report)

    cm_df = pd.DataFrame(cm, index=["Actual 0", "Actual 1"], columns=["Predicted 0", "Predicted 1"])
    logging.info("Confusion Matrix :\n%s", cm_df.to_string())

    logging.info("Model Predictions vs Actual Labels (first 50 rows):")
    for i in range(min(50, len(y_pred))):
        logging.info("Sample %d - Predicted: %s, Actual: %s", i+1, y_pred[i], y_test.iloc[i])

    print("Evaluation complete. Results saved to results.log")