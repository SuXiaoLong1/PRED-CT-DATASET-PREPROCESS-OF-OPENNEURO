# OPENNEURO的EEG数据集处理

由于OPENNEURO的EEG数据集不同数据集有不同参数。请在运行每个程序前设置参数。

下载数据集后，调整：

```

```

## 1.chul.m

chul.m将不规则的set文件转换为.mat文件

```
sub-001.mat

sub-002.mat

.....

sub-122.mat
```

## 2.chuli.m

chuli.m对文件进行分段，另存为Standard_input_nontarget_1_122.mat

```
data_path = 'F:\python\EEG\data\ds003474-download';%sub-001.mat路径
svae_path = 'F:\python\EEG\data\ds003474-download';%Standard_input_nontarget_1_122.mat路径
```

可以通过调整epoch_length调整分块长度:

```
data_length=epoch_length*fs_down
```

## 3.preprocess.m

preprocess.m将Standard_input_nontarget_1_122.mat的数据进行如下处理

```
1.基线矫正
2.滤波
3.多余通道剔除
4.异常通道替换
5.试次剔除
```

另存为 Proprocess_target_1_122.mat

space_time.m、feeatures.m的输入都是Proprocess_target_1_122.mat

## 4.space_time.m、feeatures.m

space_time.m、feeatures.m的输入都是Proprocess_target_1_122.mat

可以提取时域、频域特征。
