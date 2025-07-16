from sklearn.model_selection import train_test_split
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense
from tensorflow.keras.callbacks import EarlyStopping

def train_model(X, y):
    # Split the data
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42
    )

    # Define the model
    model = Sequential()
    model.add(Dense(128, activation='relu', input_shape=(X.shape[1],)))
    model.add(Dense(64, activation='relu'))
    model.add(Dense(1, activation='sigmoid'))  # Binary classification

    # Compile the model
    model.compile(optimizer='adam', 
                  loss='binary_crossentropy', 
                  metrics=['accuracy'])

    # Early stopping
    early_stop = EarlyStopping(patience=3, restore_best_weights=True)

    # Train the model
    history=model.fit(X_train, y_train,
              validation_data=(X_test, y_test),
              epochs=20,
              batch_size=256,
              callbacks=[early_stop])
    
    
    
    return model, history

    print(" Model training complete.")
