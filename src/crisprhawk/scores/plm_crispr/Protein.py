import random
import torch
import numpy as np
import torch.nn.functional as F
from torch import nn

seed = 10086
random.seed(seed)
np.random.seed(seed)
torch.manual_seed(seed)
torch.cuda.manual_seed(seed) 
torch.cuda.manual_seed_all(seed)
torch.backends.cudnn.benchmark = False
torch.backends.cudnn.deterministic = True

# device = torch.device("cuda:0" if torch.cuda.is_available() else 'cpu')
# torch.cuda.set_device(device)
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

class GlobalMaxPool1d(nn.Module):
    """global max pooling""" 
    def __init__(self):
        super(GlobalMaxPool1d, self).__init__()
    def forward(self, x):
         # x shape: (batch_size, channel, seq_len)
         # return shape: (batch_size, channel, 1)
        return F.max_pool1d(x, kernel_size=x.shape[2])

class TextCNN(nn.Module):
    """conv->relu->pool->dropout->linear->sigmoid"""
    def __init__(self, dropout_rate, embed_size, kernel_sizes, channel_nums):
        super(TextCNN, self).__init__()
        
        self.pool = GlobalMaxPool1d()
        self.dropout = nn.Dropout(dropout_rate)
        self.convs = nn.ModuleList()  
        for c, k in zip(channel_nums, kernel_sizes):
            self.convs.append(nn.Conv1d(in_channels=embed_size, 
                                        out_channels=c, 
                                        kernel_size=k))

    def forward(self, inputs):
        #print("input shape",inputs.shape)
        # inputs shape: (batch_size, seq_len, embed_size)
        encoding = torch.cat([self.pool(F.relu(conv(inputs.permute(0, 2, 1).float()))).squeeze(-1) for conv in self.convs], dim=1)
        return encoding

def load_protein_features(file_path):
   
    Protein_data = []
    loaded_tensor = torch.load(file_path)
    Protein_name = loaded_tensor['label']
    Protein_data = loaded_tensor['representations']

    print(Protein_name,Protein_data.size())
    return Protein_data

def extract_features(Protein_data, model):
    with torch.no_grad():
        features = model(Protein_data.to(device)).cpu().numpy()
    return features
'''
if __name__ == '__main__':
    # load data
    Protein_data = load_protein_features('eSpCas9(1.1).esm.pt')
    print('Data has generated.')
    # model 
    embed_size = 1280
    dropout_rate, kernel_sizes, channel_nums = 0.0, [5, 9, 13], [128, 128, 128]
    model = TextCNN(dropout_rate, embed_size, kernel_sizes, channel_nums).to(device)
    
    # extract features
    features = extract_features(Protein_data, model)
    print(f'Feature shape: {features.shape}')

    feature_save_path = 'eSpCas9(1.1).csv'
    df = pd.DataFrame(features)
    df.to_csv(feature_save_path, index=False)
    print(f'Features saved to {feature_save_path}')
    '''