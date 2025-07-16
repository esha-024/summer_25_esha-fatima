import matplotlib.pyplot as plt
import seaborn as sns

def visualize_training(history):
    # Plot Accuracy
    plt.figure(figsize=(8, 5))
    plt.plot(history.history['accuracy'], label='Train Accuracy')
    plt.plot(history.history['val_accuracy'], label='Validation Accuracy')
    plt.title('Training vs Validation Accuracy')
    plt.xlabel('Epoch')
    plt.ylabel('Accuracy')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig('accuracy_plot.png')
    print("Saved: accuracy_plot.png")

    # Plot Loss
    plt.figure(figsize=(8, 5))
    plt.plot(history.history['loss'], label='Train Loss')
    plt.plot(history.history['val_loss'], label='Validation Loss')
    plt.title('Training vs Validation Loss')
    plt.xlabel('Epoch')
    plt.ylabel('Loss')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig('loss_plot.png')
    print(" Saved: loss_plot.png")

def plot_feature_correlation(X, feature_names=None):
    try:
        import pandas as pd
        if not isinstance(X, pd.DataFrame):
            X = pd.DataFrame(X, columns=feature_names)

        corr = X.corr()

        plt.figure(figsize=(12, 10))
        sns.heatmap(corr, cmap='coolwarm', center=0, square=True, linewidths=0.5)
        plt.title("Feature Correlation Heatmap")
        plt.tight_layout()
        plt.savefig("feature_correlation.png")
        print(" Saved: feature_correlation.png")
    except Exception as e:
        print(f" Could not generate heatmap: {e}")
