import torch
import torch.nn as nn
import pandas as pd
import torch.nn.functional as F
import numpy as np
# Decoder for multi-class edge prediction
class GridDecoder(nn.Module):
    def __init__(self, hidden_dim, num_nodes, num_layers=2, num_heads=4):
        super().__init__()
        self.transformer = nn.TransformerEncoder(
            nn.TransformerEncoderLayer(d_model=hidden_dim, nhead=num_heads, batch_first=True),
            num_layers=num_layers
        )
        self.node_pair_predictor = nn.Sequential(
            nn.Linear(2 * hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, num_classes)
        )

    def forward(self, node_embeddings):
        transformed = self.transformer(node_embeddings.unsqueeze(0)).squeeze(0)
        N = transformed.size(0)
        expanded_i = transformed.unsqueeze(1).expand(N, N, -1)
        expanded_j = transformed.unsqueeze(0).expand(N, N, -1)
        pairwise = torch.cat([expanded_i, expanded_j], dim=-1)
        logits = self.node_pair_predictor(pairwise)
        return logits  # [N, N, num_classes]

# Edge-aware Transformer Encoder
class EdgeTypeMultiHeadAttention(nn.Module):
    def __init__(self, hidden_dim, num_heads):
        super().__init__()
        self.num_heads = num_heads
        self.attn_heads = nn.ModuleList([
            nn.MultiheadAttention(embed_dim=hidden_dim, num_heads=1, batch_first=True)
            for _ in range(num_heads)
        ])

    def forward(self, x, edge_type_matrix):
        outputs = []
        x = F.layer_norm(x, x.size()[1:])
        for i, attn in enumerate(self.attn_heads):
            edge_mask = (edge_type_matrix == (i + 1)).float()
            attn_mask = (1.0 - edge_mask).masked_fill(edge_mask == 0, -1e9)
            attn_output, _ = attn(x, x, x, attn_mask=attn_mask)
            attn_output = F.dropout(attn_output, p=0.05, training=self.training)
            outputs.append(attn_output)
        return torch.cat(outputs, dim=-1)

class DirectionalAttentionEncoder(nn.Module):
    def __init__(self, node_dim, hidden_dim, num_layers, num_heads, num_nodes):
        super().__init__()
        self.embedding = nn.Linear(node_dim, hidden_dim)
        self.attn_blocks = nn.ModuleList([
            nn.ModuleDict({
                "attn": EdgeTypeMultiHeadAttention(hidden_dim, num_heads),
                "proj": nn.Linear(hidden_dim * num_heads, hidden_dim)
            })
            for _ in range(num_layers)
        ])

    def forward(self, node_features, edge_type_matrix):
        x = self.embedding(node_features)
        for block in self.attn_blocks:
            x = block["attn"](x, edge_type_matrix)
            x = block["proj"](x)
        return x
# Define the hidden dimension
hidden_dim = 128

# Define the fully connected layer to project encoder output
fc_projection = nn.Sequential(
    nn.Linear(hidden_dim, hidden_dim),
    nn.ReLU(),
    nn.LayerNorm(hidden_dim)
)



def generate():
    symb_count = 22  # You can change this as needed
    original = np.empty((16, 16), dtype=object)
    for r in range(16):
        for c in range(16):
            vect = [-1, -1, -1, -1]  # [top, right, bottom, left]
            # Boundary conditions
            if r == 0:
                vect[0] = 0
            elif r == 15:
                vect[2] = 0
            if c == 0:
                vect[3] = 0
            elif c == 15:
                vect[1] = 0
            # Copy from neighbors
            if c > 0:
                ob = original[r, c - 1]
                vect[3] = ob[1]
            if r > 0:
                ob = original[r - 1, c]
                vect[0] = ob[2]
            # Replace remaining -1s with random integers
            vect = [np.random.randint(1, symb_count + 1) if x == -1 else x for x in vect]
            original[r, c] = vect
    return original

def generate_pair_mat(labels_tensor):
    pairs_mat = np.zeros((256, 256), dtype=int)
    #Define directions and corresponding side values
    directions = {
        (-1, 0): 1,  # up
        (0, 1): 2,  # right
        (1, 0): 3,  # down
        (0, -1): 4  # left
    }
    #Iterate over each element in the labels_tensor
    for r in range(16):
        for c in range(16):
            label1 = labels_tensor[r, c].item()
            for (dr, dc), side in directions.items():
                nr, nc = r + dr, c + dc
                if 0 <= nr < 16 and 0 <= nc < 16:
                    label2 = labels_tensor[nr, nc].item()
                    pairs_mat[label1, label2] = side

    return pairs_mat

def update_pairs_mat(pairs_mat, node_features_df, epoch_fraction):
    # Ensure using NumPy arrays
    pairs_mat = np.array(pairs_mat)
    # Loop over all pairs of vectors
    for label1 in range(len(node_features_df)):
        vector1 = np.array(node_features_df[label1])#.item())
        for label2 in range(len(node_features_df)):
            if label1 == label2:
                continue
            vector2 = np.array(node_features_df[label2])#.item())
            # Generate a shuffled testvector of [0, 1, 2, 3]
            testvector = np.random.permutation(4)
            # Loop through testvector
            for i in testvector:
                if vector1[i] == vector2[(i + 2) % 4] and pairs_mat[label1, label2] == 0 and np.random.rand()<epoch_fraction:
                    pairs_mat[label1, label2] = i + 1
                    pairs_mat[label2, label1] = (i+2)%4 + 1 #optionally: make the matrix symmetric so not to give a clue for the NN on the actual connections in the final solution
                    break
    return pairs_mat


# Parameters
grid_dim = 16
num_nodes = 256  # Input graph has exactly 256 nodes
node_dim = 4
hidden_dim = 128
num_layers = 4   #embed_dim must be divisible by num_heads
num_heads = 4
num_classes = 5  # 0: no edge, 1–4: edge types

# Instantiate models
encoder = DirectionalAttentionEncoder(node_dim, hidden_dim, num_layers, num_heads, num_nodes)
decoder = GridDecoder(hidden_dim, num_nodes, num_layers, num_heads=num_heads)

# CUDA device setup
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
encoder.to(device)
decoder.to(device)
#class_weights = class_weights.to(device)

#optimizer = torch.optim.Adam(list(encoder.parameters()) + list(decoder.parameters()), lr=3e-5)
optimizer = torch.optim.Adam(
    list(encoder.parameters()) + list(decoder.parameters()) + list(fc_projection.parameters()),
    lr=3e-5
)

#criterion = nn.CrossEntropyLoss(weight=class_weights)
#criterion = nn.BCEWithLogitsLoss()
encoder.load_state_dict(torch.load(r"C:\Users\seifer\Documents\Matlab\eternity\encoder_graph2grid.pth"))
decoder.load_state_dict(torch.load(r"C:\Users\seifer\Documents\Matlab\eternity\decoder_graph2grid.pth"))
fc_projection.load_state_dict(torch.load(r"C:\Users\seifer\Documents\Matlab\eternity\fc_projection_graph2grid.pth"))
encoder.eval()
decoder.eval()
fc_projection.eval()

# Generate synthetic graph
use_files=True
if (use_files):
    # Load node features from CSV
    node_features_df = pd.read_csv(r"C:\Users\seifer\Documents\Matlab\eternity\node_features.csv", header=None)
    node_features = torch.tensor(node_features_df.values, dtype=torch.float32).to(device)  # shape [256, node_dim]
    # Load edge type matrix from CSV
    edge_type_matrix_df = pd.read_csv(r"C:\Users\seifer\Documents\Matlab\eternity\edge_type_matrix.csv", header=None)
    edge_type_matrix = torch.tensor(edge_type_matrix_df.values, dtype=torch.float).to(device)  # shape [256, 256]
else:
    board = generate()  # temporarily fix the board
    flattened_board = board.reshape((256, 1))
    labels = np.arange(0, 256)  #between 0 and 255
    np.random.shuffle(labels)  #Now labels is shuffled
    sorted_indices = np.argsort(labels)
    sorted_vectors = flattened_board[sorted_indices]
    node_features = np.array([v[0] if isinstance(v, (list, np.ndarray)) else v for v in sorted_vectors])
    labels_reshaped = labels.reshape((16, 16))
    pairs_mat = generate_pair_mat(labels_reshaped)
    updated_pairs_mat = update_pairs_mat(pairs_mat, node_features, 0.005)
    node_features=torch.tensor(node_features, dtype=torch.float32).to(device)
    edge_type_matrix=torch.tensor(updated_pairs_mat, dtype=torch.float).to(device)
    initial_edge_type_matrix=torch.tensor(pairs_mat, dtype=torch.float).to(device)

#threshold=0.01
# Encode and decode
with torch.no_grad():
    dummy_features = torch.zeros((num_nodes, 4)).to(device)
    memory = encoder(node_features, edge_type_matrix)  # here forward is used,  consider dummy_features
    memory = fc_projection(memory)
    logits = decoder(memory)  # [256, 256, 5]
    edge_logits = logits[..., 1:].sum(dim=-1)  # Aggregate logits for edge classes
    # torch.sigmoid(edge_logits) > 0.5

    #expected_edge_count = max(960,round((edge_type_matrix > 0).sum().item()*0.2))  # or whatever you expect
    expected_edge_count = 961  # or whatever you expect
    flat_logits = edge_logits.view(-1)
    flat_edge_type_matrix=edge_type_matrix.reshape(-1)
    #threshold = torch.topk(flat_logits, expected_edge_count).values.min()
    threshold = torch.topk(flat_logits*(flat_edge_type_matrix>0), expected_edge_count).values.min()
    #predicted_clean_edge_type_matrix = (edge_logits > threshold).int()
    predicted_clean_edge_type_matrix_yesno = (edge_logits > threshold).int()
    clean_edge_type_matrix=(edge_type_matrix.int()*predicted_clean_edge_type_matrix_yesno).cpu().numpy()
    #clean_edge_type_matrix=(edge_type_matrix.int()).cpu().numpy()

# Print the predicted grid
#print(f"Predicted clean pairing matrix: {(predicted_clean_edge_type_matrix>0).sum():.4f}  Nonzero elements")
if (use_files):
    print(f"Predicted edge count: {(predicted_clean_edge_type_matrix_yesno > 0).sum().item()},  actual terms: {(edge_type_matrix*predicted_clean_edge_type_matrix_yesno > 0).sum().item()} ")
else:
    print(f"Predicted edge count: {(predicted_clean_edge_type_matrix_yesno > 0).sum().item()}, correct predictions:  {((predicted_clean_edge_type_matrix_yesno>0) & (initial_edge_type_matrix>0)).sum().item()} ")

# Save to CSV
np.savetxt(r"C:\Users\seifer\Documents\Matlab\eternity\clean_edge_type_matrix.csv", clean_edge_type_matrix, fmt='%d', delimiter=",")

