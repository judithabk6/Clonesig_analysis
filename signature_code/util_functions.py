#!/usr/bin/env python
# -*- coding:utf-8 -*-

from __future__ import absolute_import
import errno
import os


def safe_mkdir(path):
    """
    this function takes a path as input, and creates a directory without raising
    if the directory exists (and displaying that it exists)
    """
    try:
        os.mkdir(path)
    # except FileExistsError:  # for python 3.3 and more
    except OSError as e:
        if e.errno == errno.EEXIST:
            print('Directory {} not created. already exists'.format(path))
        else:
            raise
