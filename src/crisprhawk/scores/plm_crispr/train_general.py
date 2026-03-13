import torch
import torch.nn as nn
import torch.nn.functional as F

from .Protein import TextCNN

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
# device=0
# torch.cuda.set_device(device)


class sgrna_net(nn.Module):
    def __init__(self, type=11):
        super().__init__()
        self.conv1 = nn.Conv1d(
            in_channels=5, out_channels=64, kernel_size=5, padding=2
        )  # 第二条路
        self.conv2 = nn.Conv1d(
            in_channels=64, out_channels=64, kernel_size=5, padding=2
        )  # 第二条路
        self.conv3 = nn.Conv1d(
            in_channels=64, out_channels=64, kernel_size=5, padding=2
        )  # 第二条路
        self.conv4 = nn.Conv1d(
            in_channels=5, out_channels=64, kernel_size=5, padding=2
        )  # 第一条路：一层卷积

        # TextCNN for processing protein sequences
        embed_size = 1280
        self.text_cnn = TextCNN(
            dropout_rate=0.0,
            embed_size=embed_size,
            kernel_sizes=[5, 9, 13],
            channel_nums=[128, 128, 128],
        )
        # self.text_cnn = TextCNN(dropout_rate, embed_size, kernel_sizes, channel_nums)

        self.lin1 = nn.Linear(64 * 59 * 2, 2048)
        # self.lin2=nn.Linear(2048,2048)
        self.lin3 = nn.Linear(2048, 256)
        self.lin4 = nn.Linear(384, 256)  # 蛋白质处理
        self.lin5 = nn.Linear(256, 128)  # protein+sgrna
        self.lin6 = nn.Linear(128, 1)  #
        # self.bn=nn.BatchNorm1d(2048)
        # self.bn2=nn.BatchNorm1d(16)
        self.weight_net = nn.Sequential(
            nn.Linear(256 * 2, 128),
            nn.ReLU(),
            nn.Linear(128, 2),
            nn.Softmax(dim=1),  # 保证权重的和为1
        )

    def forward(self, seq, typ, protein, s1=False, train=True):
        # 蛋白质处理
        seq_protein = protein  # [N,1368,1280]
        # print("1protein shape",seq_protein.shape)
        feature = self.text_cnn(
            seq_protein.to(device)
        )  # textcnn(batch_size, channel, seq_len)
        # print("2feature",feature.shape)#(N,384）
        # print("Feature size after TextCNN:", feature.size())  # 打印特征尺寸
        # out = F.max_pool1d(feature, kernel_size=feature.shape[2])#GlobalMaxPool1d
        out_protein = F.relu(self.lin4(feature.float()))  # 蛋白质(7070,256)
        # print("3protein",out_protein.shape)
        out_protein = torch.sigmoid(out_protein)
        # print("4protein",out_protein.shape)
        if train:
            out_protein = F.dropout(out_protein, p=0.5)
        # out=self.embed1(seq_sgRNA) #(7070, 59, 16)
        # out=out.transpose(-1,-2)#(7070, 16, 59)
        # sgRNA处理
        # seq = torch.tensor(seq)
        seq_sgRNA = torch.nn.functional.one_hot(seq, num_classes=5).float()  # (N,59,5)
        # print("sgrna_seq:",seq_sgRNA.shape)#([128, 59, 5])
        out = seq_sgRNA.view(-1, 5, 59)  # 转置以匹配卷积层的期望输入(N, 5, 59)
        out = out.float()
        be = self.conv4(out)  # (7070, 64, 59)第一条路：一层卷积
        si = self.conv1(out)  # (7070, 64, 59)第二条路：三层卷积1
        # print(be.shape)
        # out=self.bn(out)
        out = torch.sigmoid(si)
        if train:
            out = F.dropout(out, p=0.5)
        out = self.conv2(out)  # 第二条路：三层卷积2
        # out=self.bn(out)
        out = torch.sigmoid(out)
        if train:
            out = F.dropout(out, p=0.5)
        out = self.conv3(out)  # 第二条路：三层卷积3
        # out=self.bn2(out)
        out = torch.sigmoid(out)
        if train:
            out = F.dropout(out, p=0.5)
        out = out.view(
            out.size(0), out.size(1) * out.size(2)
        )  # 第二条路(7070, 64 * 59)展平
        be4 = be.view(
            be.size(0), be.size(1) * be.size(2)
        )  # 第一条路#(7070, 64 * 59)展平
        out = torch.cat((out, be4), dim=-1)  # (7070, 64 * 2 * 59)
        out = self.lin1(out)  # (7070, 2048)
        # out=self.bn(out)
        out = torch.sigmoid(out)
        if train:
            out = F.dropout(out, p=0.5)
        # out=self.lin2(out)#(7070, 2048)
        # out=self.bn(out)
        # out=torch.sigmoid(out)
        # if train:
        #    out=F.dropout(out,p=0.5)
        out_sgrna = self.lin3(out)  # (7070, 256)
        # print("Shape of out_protein before concat:", out_protein.shape)
        # print("Shape of out_sgrna before concat:", out_sgrna.shape)
        # out_cat=torch.cat((out_sgrna, out_protein), dim=1)#蛋白质(7070,256)和sgrna(7070, 256)->(7070,512)

        # 将两个特征拼接起来以计算动态权重
        combined_features = torch.cat((out_sgrna, out_protein), dim=1)  # (7070, 512)
        weights = self.weight_net(combined_features)  # (7070, 2)

        # 拆分权重
        weight_sgrna = weights[:, 0].unsqueeze(1)  # (7070, 1)
        weight_protein = weights[:, 1].unsqueeze(1)  # (7070, 1)
        out_cat = weight_sgrna * out_sgrna + weight_protein * out_protein  # (7070, 256)
        out_cat = self.lin5(out_cat)  # (7070, 128)
        out_cat = torch.sigmoid(out_cat)
        if train:
            out_cat = F.dropout(out_cat, p=0.5)
        out_cat = self.lin6(out_cat)  # (7070, 1)
        out = torch.sigmoid(out_cat)  # (7070,)
        # print("Shape of output:", out.shape)
        return out
