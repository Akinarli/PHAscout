import csv

y_true = []
y_pred = []

with open('benchmark_results.csv', 'r', encoding='utf-8-sig') as f:
    reader = csv.DictReader(f)
    for row in reader:
        # Binary classification: Producer (1) vs Non-Producer (0)
        expected = 0 if row['Expected_Type'] == 'None' else 1
        predicted = 0 if (row['Predicted_Type'] == 'None' or row['Predicted_Type'] == 'belirsiz') else 1
        y_true.append(expected)
        y_pred.append(predicted)

TP = sum(1 for yt, yp in zip(y_true, y_pred) if yt == 1 and yp == 1)
TN = sum(1 for yt, yp in zip(y_true, y_pred) if yt == 0 and yp == 0)
FP = sum(1 for yt, yp in zip(y_true, y_pred) if yt == 0 and yp == 1)
FN = sum(1 for yt, yp in zip(y_true, y_pred) if yt == 1 and yp == 0)

accuracy = (TP + TN) / len(y_true)
precision = TP / (TP + FP) if (TP + FP) > 0 else 0
recall = TP / (TP + FN) if (TP + FN) > 0 else 0
f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0

print(f"Total: {len(y_true)}")
print(f"TP: {TP}, TN: {TN}, FP: {FP}, FN: {FN}")
print(f"Accuracy: {accuracy:.4f}")
print(f"Precision: {precision:.4f}")
print(f"Recall: {recall:.4f}")
print(f"F1-Score: {f1:.4f}")
