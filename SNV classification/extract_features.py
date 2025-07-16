import pandas as pd
from sklearn.preprocessing import OneHotEncoder

def extract_features(input_csv='clinvar_parsed.csv', max_total_columns=120):
    df = pd.read_csv(input_csv, low_memory=False)

    # Limit cardinality smartly
    categorical_cols = ['review_status', 'variant_class', 'consequence', 'gene']
    top_n_map = {}
    remaining_columns = max_total_columns

    for col in categorical_cols:
        unique_vals = df[col].nunique()
        allow = min(remaining_columns, unique_vals, 40)
        top_vals = df[col].value_counts().nlargest(allow).index
        df[col] = df[col].where(df[col].isin(top_vals), other='other')
        top_n_map[col] = allow
        remaining_columns -= allow
        if remaining_columns <= 0:
            break

    print(f" Encoding limits per column: {top_n_map}")

    # One-hot encode
    encoder = OneHotEncoder(sparse_output=False, handle_unknown='ignore')
    encoded = encoder.fit_transform(df[categorical_cols])
    encoded_df = pd.DataFrame(encoded, columns=encoder.get_feature_names_out(categorical_cols))

    # Numeric column
    numeric_df = df[['af_exac']].fillna(0)

    # Combine features
    X = pd.concat([encoded_df, numeric_df], axis=1)
    y = df['label']

    print(f"Extracted features with shape: {X.shape}")
    return X, y
