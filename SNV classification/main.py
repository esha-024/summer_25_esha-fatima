from extract_features import extract_features
from model import train_model
from visualize import visualize_training, plot_feature_correlation

def run_pipeline():
    print(" Extracting features...")
    X, y = extract_features()
    feature_names = list(X.columns)
    print(" Training model...")
    model, history = train_model(X, y)

    print(" Visualizing results...")
    visualize_training(history)
    plot_feature_correlation(X, feature_names)

    print(" Pipeline complete.")

if __name__ == "__main__":
    run_pipeline()
