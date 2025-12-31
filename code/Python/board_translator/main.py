# Convert Eternity II partial solution to editor format

import os
os.chdir("C:/Users/seifer/Documents/httpswww.shortestpath.seeiiresults.html")
# Convert Eternity II partial solution using e2pieces.txt (N E S W format, piece number = line number)
input_file = "467d.txt"
output_file = input_file.replace(".txt", "_converted.txt")

lookuptable = [0, 1, 6, 11, 16, 4, 9, 14, 19, 2, 7, 12, 17, 5, 10, 15, 20, 3, 8, 13, 18, 21, 22]
# Step 1: Read e2pieces.txt and build mapping of piece_number -> [N, E, S, W]
piece_map = {}
with open("e2pieces.txt", "r") as f:
    for idx, line in enumerate(f, start=1):  # start=1 for piece numbers
        parts = line.strip().split()
        if not parts:
            continue
        # Fill missing values with 0
        edges = [int(parts[i]) if i < len(parts) else 0 for i in range(0, 4)]
        piece_map[idx] = [lookuptable[edges[0]], lookuptable[edges[3]], lookuptable[edges[1]], lookuptable[edges[2]]]

# Step 2: Read 467a.txt and parse grid of piece_number/orientation
with open(input_file, "r") as f:
    text = f.read()
# Replace all occurrences of ---/- with 0/0
text = text.replace("---/-", "0/0")
# Split into lines and strip empty ones
grid_lines = [line.strip() for line in text.splitlines() if line.strip()]

grid = []
for line in grid_lines:
    row = []
    for cell in line.split():
        piece, orient = cell.split("/")
        row.append((int(piece), int(orient)))
    grid.append(row)

# Step 3: Rotation function
# Orientation: 0 -> [N,E,S,W], 1 -> [W,N,E,S], 2 -> [S,W,N,E], 3 -> [E,S,W,N]
def rotate(edges, orientation):
    if orientation == 0:
        return edges
    elif orientation == 1:
        return [edges[3], edges[0], edges[1], edges[2]]
    elif orientation == 2:
        return [edges[2], edges[3], edges[0], edges[1]]
    elif orientation == 3:
        return [edges[1], edges[2], edges[3], edges[0]]

# Step 4: Convert grid (reverse rows if needed)
converted_lines = []
for row in grid:
    converted_row = []
    for piece, orient in row:
        if piece == 0:  # Empty slot
            converted_row.append("0 0 0 0")
        else:
            rotated_edges = rotate(piece_map[piece], orient)
            converted_row.append(" ".join(map(str, rotated_edges)))
    converted_lines.append(" ".join(converted_row))  # reverse row for editor alignment

# Step 5: Save output
with open(output_file, "w") as f:
    f.write("\n".join(converted_lines))

print("Conversion complete. Output saved to converted_solution.txt")