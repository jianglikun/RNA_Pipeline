#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# @PROJECT : rna-seq-pipeline
# @Time    : 2018/3/30 15:08
# @Author  : Chen Yuelong
# @Mail    : yuelong.chen@oumeng.com.cn
# @File    : generateMakefile.py
# @Software: PyCharm

from __future__ import absolute_import, unicode_literals
import sys, os
import argparse
from src_py.utils import generate_whole_makefils

def get_args():
    '''
    命令行参数
    :return:args
    '''
    parser = argparse.ArgumentParser(description='***RNA-seq analysis pipeline(makefile) GENERATE***\n')
    parser.add_argument('--project', '-p', dest='project', required=True,
                        help='分析样本所属项目名称，详情参见SOP')
    parser.add_argument('--input_list','-l',dest='input_list',required=True,
                        help='分析样本list，详情参见SOP')
    parser.add_argument('--outdir', '-o', dest='outdir', required=True,
                        help='分析结果输出路径，详情参见SOP')
    parser.add_argument('--makefile', '-m', dest='makefile', required=True,
                        help='makefile生成路径，详情参见SOP')
    # parser.add_argument('--config', '-c', dest='config', required=True,
    #                     help='makefile所需要的config，详情参见SOP')
    args = parser.parse_args()
    return args

sys.path.append(os.path.abspath(os.path.dirname(os.path.dirname(__file__))))


def main():
    '''
    测试流程
    '''
    generate_whole_makefils(get_args())


if __name__ == '__main__':
    main()