# starrail-compass

解谜《崩坏：星穹铁道》引航罗盘的脚本。

A compass puzzle solver in HonKai: Star Rail.

## 环境

- Python 2.7.18

## 使用

需要改的是最后 `if __name__ == "__main__":` 的 `matrix` 矩阵。

前三列分别为三种转动方案，最后一列 `mod` 减号后的数字为当前刻度。
