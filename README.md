# Spectral Clustering on Multi-scale Data
This is an implementation of CAST in the <b>KDD2020</b> paper:

`CAST: A Correlation-based Adaptive Spectral Clustering Algorithm on Multi-scale Data`

> by Xiang Li, Ben Kao, Caihua Shan, Dawei Yin, Martin Ester.

Also, the project includes the implementation of ROSC in the <b>WWW2018</b> paper:

`ROSC: robust spectral clustering on multi-scale data`

> by Xiang Li, Ben Kao, Siqiang Luo, Martin Ester.

Note that both ROSC and CAST study robust spectral clustering on multi-scale data.

## Introduction

The codes are written by Matlab

The entry file is `main.m`

`ROSC.m` implements the core of ROSC

`ROSC-S.m` implements the core of ROSC-S

`CAST.m` implements the core of CAST

## Reference

Please cite the papers if you use our codes in this repo.

```
@inproceedings{li2018rosc,
  title={ROSC: Robust spectral clustering on multi-scale data},
  author={Li, Xiang and Kao, Ben and Luo, Siqiang and Ester, Martin},
  booktitle={WWW},
  pages={157--166},
  year={2018}
}

@inproceedings{li2020cast,
  title={CAST: A Correlation-based Adaptive Spectral Clustering Algorithm on Multi-scale Data},
  author={Li, Xiang and Kao, Ben and Shan, Caihua and Yin, Dawei and Ester, Martin},
  booktitle={KDD},
  year={2020}
}
```


