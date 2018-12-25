#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# @PROJECT : rna-seq-pipeline
# @Time    : 2018/3/29 10:06
# @Author  : Chen Yuelong
# @Mail    : yuelong.chen@oumeng.com.cn
# @File    : reduceDim.py
# @Software: PyCharm


'''
基因表达相关研究中经常会用到降维，所以准备将降维的方法整理一下，写到这个里面
'''


from __future__ import absolute_import, unicode_literals
import sys, os
import pandas as pd
from sklearn import manifold
sys.path.append(os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

class matrix():
    '''
    多维表达数据类
    '''
    def __init__(self,DataFrame):
        '''
        构造函数
        :param DataFrame:pandas.dataframe
        '''


def main():
    '''
    测试流程
    '''
    pass


if __name__ == '__main__':
    main()