import pandas as pd
import matplotlib.pyplot as plt
import sys

def main():
    if len(sys.argv) < 2:
        print("Usage: python plot_symmetry.py output_filename.png")
        sys.exit(1)

    output_filename = sys.argv[1]

    # Read the fixed CSV file
    df = pd.read_csv("implicit_symmetry_focal_set.csv")
    curve = df[df['type'] == 'curve']
    symmetry = df[df['type'] == 'symmetry']
    focal = df[df['type'] == 'focal']

    # Plot
    plt.figure(figsize=(10, 10))
    plt.plot(curve['x'], curve['y'], 'b-', linewidth=2, label='Curve')
    plt.scatter(symmetry['x'], symmetry['y'], s=3, color='red', alpha=0.6, label='Symmetry Set')
    plt.scatter(focal['x'], focal['y'], s=2, color='green', alpha=0.6, label='Focal Set')
    plt.axis('equal')
    plt.legend()
    plt.title("Curve and Symmetry Set")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.grid(True)

    # Save and show
    plt.savefig(output_filename, dpi=300)
    print(f"Saved plot to {output_filename}")
    plt.show()

if __name__ == "__main__":
    main()
