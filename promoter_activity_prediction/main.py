from features_extraction import load_and_prepare, extract_all_features
from model import train_test_model
from evaluation import evaluate_and_log

def main():
    df = load_and_prepare("tableData.csv")
    X, y, df_feat = extract_all_features(df, k=4)
    model, X_test, y_test = train_test_model(X, y)
    evaluate_and_log(model, X_test, y_test, df_feat)

if __name__ == "__main__":
    main()