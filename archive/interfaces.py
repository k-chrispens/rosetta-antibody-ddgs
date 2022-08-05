"""Finding interface residues and appending a 
column to the data telling whether interface 
residues are in the mutations.
Author: Karson Chrispens"""

import re
import numpy as np
import pandas as pd
data = pd.read_csv("./raw_datasets/full_data.csv")
data.dropna(subset="Mutations", inplace=True)
# This is terribly messy but will have to do -
# used pymol script to extend function, pymol_interface.py, then calculated interface lists for all these. 
# Used regex search and replace to remove last letter of positions like 100A, 100B, 100C, etc.

pdb_1BJ1 = ['V2', 'G26', 'Y27', 'T28', 'F29', 'T30', 'N31', 'Y32', 'G33', 'M34', 'N35', 'W36', 'W47', 'V48', 'G49', 'W50', 'I51', 'N52', 'T53', 'Y54', 'T55', 'G56', 'E57', 'P58', 'T59', 'Y60', 'A61', 'A62', 'L72', 'T74', 'S77',
            'A79', 'A97', 'K98', 'Y99', 'P100', 'H101', 'Y102', 'Y103', 'G104', 'S105', 'S106', 'H107', 'W108', 'Y109', 'F110', 'D111', 'V112', 'D1', 'I29', 'Y32', 'Y49', 'F50', 'Q90', 'Y91', 'S92', 'T93', 'V94', 'P95', 'W96', 'T97']
pdb_1NMB = ['T28', 'F29', 'T30', 'N31', 'Y32', 'N33', 'M34', 'Y35', 'G49', 'I50', 'F51', 'Y52', 'P52', 'G53', 'N54', 'G55', 'D56', 'T57', 'S58', 'Y59', 'N60', 'Q61', 'T68', 'L69', 'T70', 'A71', 'D72', 'K73', 'S74', 'S95', 'G96', 'G97', 'S98', 'Y99',
            'R100', 'Y100', 'D100', 'G100', 'G100', 'F100', 'D1', 'I2', 'Q3', 'A25', 'S26', 'Q27', 'D28', 'I29', 'S30', 'N31', 'Y32', 'L33', 'N34', 'Y49', 'Y50', 'T51', 'G66', 'S67', 'G68', 'T69', 'Q89', 'Q90', 'D91', 'F92', 'T93', 'L94', 'P95', 'F96', 'T97']
pdb_3BDY = ['K30', 'D31', 'T32', 'Y33', 'I34', 'H35', 'A49', 'R50', 'I51', 'Y52', 'P52', 'T53', 'N54', 'G55', 'Y56', 'T57', 'R58', 'Y59', 'R94', 'W95', 'G96', 'G97', 'D98', 'G99', 'F100', 'Y100', 'A101', 'M102', 'D103', 'D1', 'I2', 'A25', 'Q27',
            'I27', 'P27', 'R27', 'S28', 'I29', 'S30', 'G31', 'Y32', 'V33', 'A34', 'I48', 'Y49', 'W50', 'G51', 'S52', 'Y53', 'L54', 'G64', 'S65', 'G66', 'S67', 'G68', 'T69', 'D70', 'F71', 'Q89', 'Q90', 'H91', 'Y92', 'T93', 'T94', 'P95', 'P96', 'T97']
pdb_1C08 = ['D1', 'I2', 'Q27', 'S28', 'I29', 'G30', 'N31', 'N32', 'L33', 'H34', 'I48', 'K49', 'Y50', 'A51', 'S52', 'Q53', 'S54', 'I55', 'S56', 'S65', 'G66', 'S67', 'G68', 'T69', 'Q89', 'Q90', 'S91', 'N92', 'S93', 'W94', 'P95', 'Y96', 'T97', 'D1',
            'V2', 'G26', 'D27', 'S28', 'I29', 'T30', 'S31', 'D32', 'Y33', 'W34', 'S35', 'Y47', 'M48', 'G49', 'Y50', 'V51', 'S52', 'Y53', 'S54', 'G55', 'S56', 'T57', 'Y58', 'Y59', 'N60', 'P61', 'T73', 'S74', 'A96', 'N97', 'W98', 'D99', 'G100', 'D101', 'Y102']
pdb_1VFB = ['Q1', 'V2', 'G26', 'F27', 'S28', 'L29', 'T30', 'G31', 'Y32', 'G33', 'V34', 'N35', 'M50', 'I51', 'W52', 'G53', 'D54', 'G55', 'N56', 'T57', 'D58', 'R97', 'E98', 'R99', 'D100', 'Y101', 'R102', 'L103', 'D104', 'D1', 'I2',
            'G27', 'N28', 'I29', 'H30', 'N31', 'Y32', 'L33', 'A34', 'L47', 'V48', 'Y49', 'Y50', 'T51', 'T52', 'T53', 'L54', 'A55', 'D56', 'G64', 'S65', 'G66', 'S67', 'G68', 'Q89', 'H90', 'F91', 'W92', 'S93', 'T94', 'P95', 'R96', 'T97']
pdb_3BE1 = ['G16', 'S17', 'K30', 'D31', 'T32', 'Y33', 'I34', 'H35', 'A49', 'R50', 'I51', 'Y52', 'P52', 'T53', 'N54', 'G55', 'Y56', 'T57', 'R58', 'Y59', 'A60', 'D61', 'K64', 'G65', 'N82', 'R94', 'W95', 'G96', 'G97', 'D98', 'G99', 'F100', 'Y100', 'A101', 'M102',
            'D103', 'D1', 'I2', 'D27', 'I27', 'P27', 'R27', 'S28', 'I29', 'S30', 'G31', 'Y32', 'V33', 'A34', 'I48', 'Y49', 'W50', 'G51', 'S52', 'Y53', 'L54', 'Y55', 'S56', 'G66', 'S67', 'G68', 'T69', 'Q89', 'Q90', 'H91', 'Y92', 'T93', 'T94', 'P95', 'P96', 'T97']
pdb_1CZ8 = ['V2', 'G26', 'Y27', 'D28', 'F29', 'T30', 'H31', 'Y32', 'G33', 'M34', 'N35', 'W36', 'W47', 'V48', 'G49', 'W50', 'I51', 'N52', 'T53', 'Y54', 'T55', 'G56', 'E57', 'P58',
            'T59', 'Y60', 'A61', 'A62', 'L72', 'T74', 'A79', 'A97', 'K98', 'Y99', 'P100', 'Y101', 'Y102', 'Y103', 'G104', 'T105', 'S106', 'H107', 'W108', 'Y109', 'F110', 'D111', 'V112']
pdb_1XGP = ['D1', 'I2', 'Q27', 'S28', 'I29', 'S30', 'N31', 'N32', 'L33', 'H34', 'I48', 'K49', 'Y50', 'A51', 'S52', 'Q53', 'S54', 'I55', 'S56', 'G57', 'S65', 'G66', 'S67', 'G68', 'T69', 'Q89', 'Q90', 'S91', 'N92', 'S93', 'W94', 'P95', 'Y96', 'T97', 'E1',
            'V2', 'Q3', 'G26', 'D27', 'S28', 'V29', 'T30', 'S31', 'D32', 'A33', 'W34', 'S35', 'Y47', 'M48', 'G49', 'Y50', 'I51', 'S52', 'Y53', 'S54', 'G55', 'S56', 'T57', 'Y58', 'Y59', 'H60', 'P61', 'T73', 'S74', 'A96', 'S97', 'W98', 'G99', 'G100', 'D101', 'V102']
pdb_3BN9 = ['D1', 'I2', 'A25', 'Q27', 'G28', 'I29', 'S30', 'S31', 'Y32', 'L33', 'Y49', 'A50', 'A51', 'S52', 'S53', 'S67', 'G68', 'T69', 'Q90', 'H91', 'G92', 'N93', 'L94', 'P95', 'Y96', 'T97', 'S25', 'G26', 'F27', 'T28', 'F29', 'S30', 'S31', 'Y32', 'A33', 'M34', 'S35', 'W36', 'W47', 'V48', 'S49',
            'A50', 'I51', 'S52', 'G52', 'S53', 'G54', 'G55', 'S56', 'T57', 'Y58', 'Y59', 'A60', 'I69', 'S70', 'R71', 'D72', 'N73', 'S74', 'K75', 'N76', 'A93', 'R94', 'P95', 'Y96', 'L97', 'T98', 'Y99', 'P100', 'Q100', 'R100', 'R100', 'G100', 'P100', 'Q100', 'N100', 'V100', 'S100', 'P100', 'F100']
pdb_1DQJ = ['E1', 'V2', 'Q3', 'G26', 'D27', 'S28', 'V29', 'T30', 'S31', 'D32', 'Y33', 'W34', 'S35', 'Y47', 'M48', 'G49', 'Y50', 'I51', 'S52', 'Y53', 'S54', 'G55', 'S56', 'T57', 'Y58', 'Y59', 'H60', 'P61', 'T73', 'S74', 'N76', 'A96', 'S97', 'W98', 'G99',
            'G100', 'D101', 'V102', 'D1', 'I2', 'Q27', 'S28', 'I29', 'S30', 'N31', 'N32', 'L33', 'H34', 'I48', 'K49', 'Y50', 'A51', 'S52', 'Q53', 'S54', 'I55', 'S56', 'G57', 'S65', 'G66', 'S67', 'G68', 'T69', 'Q89', 'Q90', 'S91', 'N92', 'S93', 'W94', 'P95', 'Y96', 'T97']
pdb_1XGQ = ['D1', 'I2', 'Q27', 'S28', 'I29', 'S30', 'N31', 'N32', 'L33', 'H34', 'I48', 'K49', 'Y50', 'A51', 'S52', 'Q53', 'S54', 'I55', 'S56', 'S65', 'G66', 'S67', 'G68', 'T69', 'Q89', 'Q90', 'S91', 'N92', 'S93', 'W94', 'P95', 'Y96', 'T97', 'E1', 'V2',
            'Q3', 'V24', 'G26', 'D27', 'S28', 'V29', 'T30', 'S31', 'D32', 'V33', 'W34', 'S35', 'Y47', 'M48', 'G49', 'Y50', 'I51', 'S52', 'Y53', 'S54', 'G55', 'S56', 'T57', 'Y58', 'Y59', 'H60', 'P61', 'T73', 'A96', 'S97', 'W98', 'G99', 'G100', 'D101', 'V102']
pdb_3G6D = ['V2', 'G26', 'F27', 'T28', 'F29', 'N30', 'S31', 'Y32', 'W33', 'I34', 'A52', 'Y53', 'D54', 'S55', 'S56', 'N57', 'R98', 'G99', 'L100', 'G101', 'A102', 'F103', 'H104', 'W105', 'D106', 'M107', 'Q108', 'P109', 'D110', 'Y111', 'S1', 'Y2', 'S23', 'G24',
            'D25', 'N26', 'I27', 'G28', 'G29', 'T30', 'F31', 'V32', 'S33', 'I47', 'Y48', 'D49', 'D50', 'N51', 'D52', 'R53', 'P54', 'G63', 'S64', 'N65', 'S66', 'G67', 'N68', 'T69', 'A70', 'T71', 'G88', 'T89', 'W90', 'D91', 'M92', 'V93', 'T94', 'N95', 'N96', 'V97']
pdb_1DVF = ['D1', 'I2', 'G27', 'N28', 'I29', 'H30', 'N31', 'Y32', 'L33', 'A34', 'L47', 'V48', 'Y49', 'Y50', 'T51', 'T52', 'T53', 'L54', 'A55', 'D56', 'G57', 'V58', 'G64', 'S67', 'G68', 'H90', 'F91', 'W92', 'S93', 'T94', 'P95', 'R96', 'Q1', 'V2', 'G26', 'F27', 'S28', 'L29', 'T30', 'G31', 'Y32', 'G33', 'V34', 'N35', 'G49', 'M50', 'I51', 'W52', 'G53', 'D54', 'G55', 'N56', 'T57', 'D58', 'Y59', 'N60', 'I69', 'N73', 'A96', 'R97', 'E98', 'R99', 'D100', 'Y101',
            'R102', 'L103', 'D104', 'D1', 'I2', 'Q27', 'D28', 'I29', 'S30', 'N31', 'Y32', 'L33', 'N34', 'I48', 'Y49', 'Y50', 'T51', 'S52', 'R53', 'L54', 'H55', 'S56', 'G57', 'G68', 'Q90', 'G91', 'N92', 'T93', 'L94', 'P95', 'W96', 'V2', 'F27', 'N28', 'I29', 'K30', 'D31', 'T32', 'H33', 'M34', 'R50', 'I51', 'D52', 'P52', 'A53', 'N54', 'G55', 'N56', 'I57', 'Q58', 'T73', 'T94', 'K95', 'V96', 'I97', 'Y98', 'Y99', 'Q100', 'G100', 'R100', 'G100', 'A100', 'M100', 'D101']
pdb_1XGR = ['D1', 'I2', 'Q27', 'S28', 'I29', 'S30', 'N31', 'N32', 'L33', 'H34', 'I48', 'K49', 'Y50', 'A51', 'S52', 'Q53', 'S54', 'I55', 'S56', 'S65', 'G66', 'S67', 'G68', 'T69', 'Q89', 'Q90', 'S91', 'N92', 'S93', 'W94', 'P95', 'Y96', 'T97', 'A1', 'V2',
            'Q3', 'G26', 'D27', 'S28', 'V29', 'T30', 'S31', 'D32', 'I33', 'W34', 'S35', 'Y47', 'M48', 'G49', 'Y50', 'I51', 'S52', 'Y53', 'S54', 'G55', 'S56', 'T57', 'Y58', 'Y59', 'H60', 'P61', 'T73', 'S74', 'A96', 'S97', 'W98', 'G99', 'G100', 'D101', 'V102']
pdb_3GBN = ['V2', 'A24', 'S25', 'G26', 'G27', 'P28', 'F29', 'R30', 'S31', 'Y32', 'A33', 'I34', 'I51', 'I52', 'P52', 'I53', 'F54', 'G55', 'T56', 'T57', 'K58',
            'D72', 'D73', 'F74', 'A75', 'G76', 'T77', 'K94', 'H95', 'M96', 'G97', 'Y98', 'Q99', 'V100', 'R100', 'E100', 'T100', 'M100', 'N30', 'D31', 'Y32', 'D50', 'A95']
pdb_1JRH = ['L29', 'T30', 'T31', 'Y32', 'G33', 'M34', 'G35', 'V35', 'G35', 'W36', 'W47', 'L48', 'A49', 'H50', 'I51', 'W52', 'W53', 'D54', 'D55', 'D56', 'K57', 'Y58', 'Y59', 'N60', 'P61', 'S62', 'L63', 'K64', 'A93', 'R94', 'R95', 'A96', 'P97', 'F98',
            'Y99', 'G100', 'N100', 'H100', 'A100', 'M100', 'S1', 'V2', 'E3', 'A25', 'S26', 'E27', 'D28', 'I29', 'Y30', 'N31', 'R32', 'L33', 'A34', 'S49', 'G50', 'A51', 'S67', 'G68', 'K69', 'D70', 'Q89', 'Q90', 'Y91', 'W92', 'S93', 'T94', 'W96', 'T97', 'F98']
pdb_1XGT = ['D1', 'I2', 'Q27', 'S28', 'I29', 'S30', 'N31', 'N32', 'L33', 'H34', 'I48', 'K49', 'Y50', 'A51', 'S52', 'Q53', 'S54', 'I55', 'S56', 'G57', 'S65', 'G66', 'S67', 'G68', 'T69', 'Q89', 'Q90', 'S91', 'N92', 'S93', 'W94', 'P95', 'Y96', 'T97', 'E1', 'V2',
            'Q3', 'G26', 'D27', 'S28', 'V29', 'T30', 'S31', 'D32', 'L33', 'W34', 'S35', 'Y47', 'M48', 'G49', 'Y50', 'I51', 'S52', 'Y53', 'S54', 'G55', 'S56', 'T57', 'Y58', 'Y59', 'H60', 'P61', 'T73', 'S74', 'N76', 'A96', 'S97', 'W98', 'G99', 'G100', 'D101', 'V102']
pdb_3HFM = ['D1', 'V2', 'G26', 'D27', 'S28', 'I29', 'T30', 'S31', 'D32', 'Y33', 'W34', 'S35', 'Y47', 'M48', 'G49', 'Y50', 'V51', 'S52', 'Y53', 'S54', 'G55', 'S56', 'T57', 'Y58', 'Y59', 'N60', 'P61', 'T73', 'S74', 'N76', 'A96', 'N97', 'W98', 'D99', 'G100',
            'D101', 'Y102', 'D1', 'I2', 'Q27', 'S28', 'I29', 'G30', 'N31', 'N32', 'L33', 'H34', 'L46', 'I48', 'K49', 'Y50', 'A51', 'S52', 'Q53', 'S54', 'I55', 'S56', 'S65', 'G66', 'S67', 'G68', 'T69', 'Q89', 'Q90', 'S91', 'N92', 'S93', 'W94', 'P95', 'Y96', 'T97']
pdb_1KIP = ['D1', 'I2', 'G27', 'N28', 'I29', 'H30', 'N31', 'Y32', 'L33', 'A34', 'L47', 'V48', 'Y49', 'Y50', 'T51', 'T52', 'T53', 'L54', 'A55', 'D56', 'S63', 'G64', 'S65', 'G66', 'S67', 'G68', 'Q89', 'H90', 'F91', 'W92', 'S93', 'T94',
            'P95', 'R96', 'T97', 'Q1', 'V2', 'G26', 'F27', 'S28', 'L29', 'T30', 'G31', 'A32', 'G33', 'V34', 'N35', 'M50', 'I51', 'W52', 'G53', 'D54', 'G55', 'N56', 'T57', 'D58', 'R97', 'E98', 'R99', 'D100', 'Y101', 'R102', 'L103', 'D104']
pdb_1XGU = ['D1', 'I2', 'Q27', 'S28', 'I29', 'S30', 'N31', 'N32', 'L33', 'H34', 'I48', 'K49', 'Y50', 'A51', 'S52', 'Q53', 'S54', 'I55', 'S56', 'G57', 'S65', 'G66', 'S67', 'G68', 'T69', 'Q89', 'Q90', 'S91', 'N92', 'S93', 'W94', 'P95', 'Y96', 'T97', 'E1', 'V2',
            'Q3', 'G26', 'D27', 'S28', 'V29', 'T30', 'S31', 'D32', 'F33', 'W34', 'S35', 'Y47', 'M48', 'G49', 'Y50', 'I51', 'S52', 'Y53', 'S54', 'G55', 'S56', 'T57', 'Y58', 'Y59', 'H60', 'P61', 'T73', 'S74', 'N76', 'A96', 'S97', 'W98', 'G99', 'G100', 'D101', 'V102']
pdb_3L5X = ['L29', 'S30', 'T31', 'Y32', 'G33', 'M34', 'G35', 'V36', 'G37', 'A51', 'H52', 'I53', 'W54', 'W55', 'D56', 'D57', 'V58', 'K59', 'R60', 'Y61', 'N62', 'P63', 'R99', 'M100', 'G101',
            'S102', 'D103', 'Y104', 'D105', 'V106', 'W107', 'K27', 'S28', 'I29', 'S30', 'K31', 'Y32', 'L33', 'Y49', 'S50', 'G51', 'S67', 'G68', 'Q90', 'H91', 'N92', 'E93', 'Y94', 'P95', 'Y96']
pdb_1KIQ = ['D1', 'I2', 'G27', 'N28', 'I29', 'H30', 'N31', 'Y32', 'L33', 'A34', 'L47', 'V48', 'Y49', 'Y50', 'T51', 'T52', 'T53', 'L54', 'A55', 'D56', 'S63', 'G64', 'S65', 'G66', 'S67', 'G68', 'Q89', 'H90', 'F91', 'W92', 'S93',
            'T94', 'P95', 'R96', 'T97', 'V2', 'G26', 'F27', 'S28', 'L29', 'T30', 'G31', 'Y32', 'G33', 'V34', 'N35', 'M50', 'I51', 'W52', 'G53', 'D54', 'G55', 'N56', 'T57', 'D58', 'R97', 'E98', 'R99', 'D100', 'F101', 'R102', 'D104']
pdb_1YY9 = ['P14', 'S15', 'F27', 'S28', 'L29', 'T30', 'N31', 'Y32', 'G33', 'V34', 'H35', 'S40', 'G42', 'K43', 'G44', 'E46', 'W47', 'L48', 'G49', 'V50', 'I51', 'W52', 'S53', 'G54', 'G55', 'N56', 'T57', 'D58', 'Y59', 'N60', 'T61', 'P62', 'F63', 'T64', 'S65', 'R66', 'S68', 'I69', 'N70', 'K71', 'D72', 'N73', 'S74', 'V78', 'F79', 'N83', 'S84',
            'L85', 'Q86', 'S87', 'N88', 'D89', 'T90', 'A91', 'R97', 'A98', 'L99', 'T100', 'Y101', 'Y102', 'D103', 'Y104', 'E105', 'F106', 'A107', 'V117', 'D1', 'I2', 'L3', 'A25', 'S26', 'Q27', 'S28', 'I29', 'G30', 'T31', 'N32', 'I33', 'H34', 'K49', 'Y50', 'A51', 'E53', 'G68', 'T69', 'Q89', 'Q90', 'N91', 'N92', 'N93', 'W94', 'P95', 'T96', 'T97', 'F98']
pdb_3N85 = ['E1', 'V2', 'Q3', 'A24', 'S25', 'G26', 'F27', 'S28', 'I29', 'W30', 'W31', 'S32', 'W33', 'I34', 'H35', 'A49', 'S50', 'I51', 'S52', 'P52', 'S53', 'S54', 'G55', 'W56', 'T57', 'S58', 'Y59', 'A60', 'D61', 'S62', 'K64', 'T68', 'I69', 'S70', 'D72', 'T73', 'S74',
            'K75', 'N76', 'A93', 'R94', 'W95', 'W96', 'S97', 'S98', 'A99', 'M100', 'D100', 'Y101', 'D1', 'I2', 'Q27', 'S28', 'V29', 'S30', 'S31', 'A32', 'V33', 'A34', 'Y49', 'S50', 'A51', 'S52', 'S53', 'Y55', 'S56', 'Q89', 'Q90', 'W91', 'W92', 'W93', 'W94', 'P95', 'S96', 'T97']
pdb_1KIR = ['D1', 'I2', 'G27', 'N28', 'I29', 'H30', 'N31', 'Y32', 'L33', 'A34', 'V48', 'Y49', 'S50', 'T51', 'T52', 'T53', 'L54', 'A55', 'D56', 'S65', 'G66', 'S67', 'G68', 'T69', 'Q89', 'H90', 'F91', 'W92', 'S93', 'T94', 'P95',
            'R96', 'T97', 'Q1', 'V2', 'G26', 'F27', 'S28', 'L29', 'T30', 'G31', 'Y32', 'G33', 'V34', 'N35', 'M50', 'I51', 'W52', 'G53', 'D54', 'G55', 'N56', 'T57', 'D58', 'N73', 'R97', 'E98', 'R99', 'D100', 'Y101', 'R102', 'L103', 'D104']
pdb_2B2X = ['F27', 'T28', 'F29', 'S30', 'R31', 'Y32', 'T33', 'M34', 'S35', 'W36', 'V37', 'E46', 'W47', 'V48', 'A49', 'V50', 'I51', 'S52', 'G53', 'G54', 'G55', 'H56', 'T57', 'Y58', 'Y59', 'L60', 'D61', 'S62', 'V63', 'E64', 'I69', 'R71', 'D72', 'N73', 'S74', 'N76', 'R97', 'G98', 'F99', 'G100',
            'D101', 'G102', 'G103', 'Y104', 'F105', 'D106', 'Q1', 'I2', 'A25', 'S26', 'S27', 'Q28', 'V29', 'N30', 'H31', 'M32', 'F33', 'I47', 'Y48', 'L49', 'T50', 'S51', 'Y52', 'L53', 'A54', 'S55', 'G63', 'S64', 'G65', 'S66', 'G67', 'T68', 'Q88', 'Q89', 'W90', 'S91', 'G92', 'N93', 'P94', 'W95', 'T96']
pdb_3NGB = ['E28', 'F29', 'I30', 'D31', 'C32', 'T33', 'L34', 'N35', 'W36', 'I37', 'R44', 'P45', 'E46', 'W47', 'M48', 'G49', 'W50', 'L51', 'K52', 'P52', 'R53', 'G54', 'G55', 'A56', 'V57', 'N58', 'Y59', 'A60', 'R61', 'P62', 'L63', 'Q64', 'G65', 'R66', 'V67', 'T68', 'M69', 'T70', 'R71', 'D72', 'V73',
            'Y74', 'S75', 'D76', 'T77', 'A78', 'F79', 'L80', 'E81', 'R94', 'G95', 'K96', 'N97', 'C98', 'D99', 'Y100', 'N100', 'W100', 'D100', 'F100', 'V3', 'L4', 'T5', 'R24', 'T25', 'S26', 'Q27', 'Y28', 'G29', 'S30', 'L33', 'A34', 'Y49', 'S50', 'G51', 'S52', 'Q89', 'Q90', 'Y91', 'E96', 'F97', 'F98', 'G99']
pdb_1MHP = ['F27', 'T28', 'F29', 'S30', 'R31', 'Y32', 'T33', 'M34', 'S35', 'W36', 'V37', 'E46', 'W47', 'V48', 'A49', 'T50', 'I51', 'S52', 'G53', 'G54', 'G55', 'H56', 'T57', 'Y58', 'Y59', 'L60', 'D61', 'S62', 'V63', 'K64', 'I69', 'R71', 'D72', 'N73', 'S74', 'N76', 'R97', 'G98',
            'F99', 'G100', 'D101', 'G102', 'G103', 'Y104', 'F105', 'I2', 'S24', 'A25', 'S27', 'S28', 'V29', 'N30', 'H31', 'M32', 'F33', 'Y48', 'L49', 'T50', 'S51', 'N52', 'S64', 'G65', 'S66', 'G67', 'T68', 'D69', 'Y70', 'Q88', 'Q89', 'W90', 'S91', 'G92', 'N93', 'P94', 'W95', 'T96']
pdb_2BDN = ['E1', 'V2', 'Q3', 'L4', 'A24', 'S25', 'G26', 'L27', 'N28', 'I29', 'K30', 'D31', 'T32', 'Y33', 'M34', 'H35', 'R50', 'I51', 'D52', 'P53', 'A54', 'N55', 'G56', 'N57', 'T58', 'K59', 'A72', 'D73', 'T74', 'N77', 'A79', 'A97', 'R98', 'G99', 'V100', 'F101', 'G102',
            'F103', 'F104', 'D105', 'Y106', 'W107', 'E27', 'D28', 'I29', 'Y30', 'N31', 'R32', 'L33', 'A34', 'L46', 'L47', 'I48', 'S49', 'G50', 'A51', 'T52', 'S53', 'L54', 'E55', 'T56', 'G57', 'V58', 'P59', 'S60', 'R61', 'G66', 'S67', 'G68', 'K69', 'Q90', 'F91', 'W92', 'S93', 'A94']
pdb_3SE8 = ['N28', 'F29', 'R30', 'D31', 'Y32', 'S33', 'I34', 'H35', 'W36', 'V37', 'R38', 'K43', 'G44', 'F45', 'E46', 'W47', 'I48', 'G49', 'W50', 'I51', 'K52', 'P52', 'L53', 'W54', 'G55', 'A56', 'V57', 'S58', 'Y59', 'A60', 'R61', 'Q62', 'L63', 'Q64', 'G65', 'R66', 'V67', 'S68', 'M69', 'T70', 'R71', 'Q72', 'L73', 'S74', 'Q75', 'D76', 'P76', 'D76', 'D76', 'P76',
            'D76', 'W76', 'G76', 'A78', 'M80', 'E81', 'R94', 'R95', 'G96', 'S97', 'C98', 'D99', 'Y100', 'C100', 'G100', 'D100', 'F100', 'P100', 'W100', 'E1', 'I2', 'V3', 'L4', 'A25', 'S26', 'Q27', 'G29', 'G30', 'N31', 'A32', 'M33', 'T34', 'W35', 'I48', 'Y49', 'D50', 'T51', 'S52', 'R53', 'S67', 'G68', 'T69', 'C88', 'Q89', 'Q90', 'F91', 'E96', 'F97', 'F98', 'G99']
pdb_1MLC = ['D1', 'I2', 'S26', 'Q27', 'S28', 'I29', 'S30', 'N31', 'N32', 'L33', 'H34', 'K49', 'Y50', 'V51', 'Q53', 'G68', 'Q89', 'Q90', 'S91', 'N92', 'S93', 'W94', 'P95', 'R96', 'T97', 'Y27', 'T28', 'F29', 'S30', 'T31', 'Y32', 'W33',
            'I34', 'E35', 'W36', 'W47', 'I48', 'G49', 'E50', 'I51', 'L52', 'P53', 'G54', 'S55', 'G56', 'S57', 'T58', 'Y59', 'Y60', 'N61', 'E62', 'F70', 'T71', 'A72', 'T74', 'A97', 'R98', 'G99', 'D100', 'G101', 'N102', 'Y103', 'G104', 'Y105']
pdb_2NY7 = ['Q1', 'V2', 'A24', 'S25', 'G26', 'Y27', 'R28', 'F29', 'S30', 'N31', 'F32', 'V33', 'I34', 'H35', 'G49', 'W50', 'I51', 'N52', 'P52', 'Y53', 'N54', 'G55', 'N56', 'K57', 'E58', 'F59', 'S60', 'A61', 'Q64', 'F69', 'T70', 'A71',
            'D72', 'T73', 'S74', 'A75', 'N76', 'T77', 'A93', 'R94', 'V95', 'G96', 'P97', 'Y98', 'S99', 'W100', 'D100', 'D100', 'S100', 'P100', 'Q100', 'D100', 'N100', 'Y100', 'Y100', 'M100', 'D101', 'V102', 'G92', 'A93', 'S94', 'S95']
pdb_4FQY = ['V2', 'S24', 'S25', 'G26', 'G27', 'T28', 'S29', 'N30', 'N31', 'Y32', 'A33', 'I34', 'G50', 'I51', 'S52', 'P52', 'I53', 'F54', 'G55', 'S56', 'T57', 'A58', 'Y59', 'S70', 'A71', 'D72', 'I73', 'F74', 'S75', 'N76',
            'T77', 'R94', 'H95', 'G96', 'N97', 'Y98', 'Y99', 'Y100', 'Y100', 'S100', 'G100', 'M100', 'G29', 'R30', 'R31', 'S32', 'V33', 'Y49', 'S50', 'N51', 'A90', 'W91', 'D92', 'D93', 'S94', 'L95', 'K95', 'G95', 'A96']
pdb_1N8Z = ['K30', 'D31', 'T32', 'Y33', 'I34', 'H35', 'W47', 'A49', 'R50', 'I51', 'Y52', 'P53', 'T54', 'N55', 'G56', 'Y57', 'T58', 'R59', 'Y60', 'A61', 'R98', 'W99', 'G100', 'G101', 'D102', 'G103', 'F104', 'Y105', 'A106', 'M107', 'D108',
            'D1', 'I2', 'Q27', 'D28', 'V29', 'N30', 'T31', 'A32', 'V33', 'A34', 'I48', 'Y49', 'S50', 'A51', 'S52', 'F53', 'L54', 'Y55', 'S56', 'G64', 'S65', 'R66', 'S67', 'G68', 'T69', 'Q89', 'Q90', 'H91', 'Y92', 'T93', 'T94', 'P95', 'P96', 'T97']
pdb_2NYY = ['F27', 'T28', 'F29', 'K30', 'Y31', 'D32', 'Y33', 'M34', 'Y35', 'A49', 'T50', 'I51', 'S52', 'D53', 'G54', 'G55', 'S56', 'Y57', 'T58', 'Y59', 'Y60', 'S61', 'D62', 'E65', 'S71', 'R72', 'D73', 'N74', 'R98', 'Y99', 'R100', 'Y101', 'D102', 'D103', 'A104', 'M105',
            'E1', 'I2', 'V3', 'A25', 'S26', 'E27', 'S28', 'V29', 'D30', 'S31', 'Y32', 'G33', 'H34', 'S35', 'F36', 'M37', 'Q38', 'Y53', 'R54', 'A55', 'S56', 'N57', 'G68', 'S69', 'G70', 'S71', 'G72', 'T73', 'D74', 'F75', 'Q93', 'Q94', 'G95', 'N96', 'E97', 'V98', 'P99', 'F100', 'T101']
pdb_5C6T = ['Q1', 'L2', 'Q3', 'L4', 'S25', 'G26', 'A27', 'S28', 'I29', 'D30', 'R31', 'S32', 'T33', 'Y34', 'Y35', 'W36', 'I53', 'Y54', 'Y55', 'N56', 'G57', 'R58', 'A98', 'T99', 'R100', 'W101', 'N102', 'Y103', 'F104', 'F105',
            'D106', 'F107', 'D108', 'Y109', 'W110', 'G111', 'P14', 'G15', 'S35', 'A44', 'P45', 'K46', 'L47', 'L48', 'I49', 'Y50', 'R51', 'N52', 'N53', 'Q54', 'R55', 'P56', 'S57', 'G58', 'V59', 'P60', 'D61', 'S64', 'G65', 'S66']
pdb_1NCA = ['Y27', 'T28', 'F29', 'T30', 'N31', 'Y32', 'G33', 'M34', 'N35', 'G49', 'W50', 'I51', 'N52', 'T52', 'N53', 'T54', 'G55', 'E56', 'P57', 'T58', 'Y59', 'G60', 'E61', 'K64', 'T73', 'R94', 'G95', 'E96', 'D97', 'N98', 'F99', 'G100', 'S100', 'L100', 'S100', 'D101', 'Y102', 'D1',
            'I2', 'Q27', 'D28', 'V29', 'S30', 'T31', 'A32', 'V33', 'V34', 'W35', 'L46', 'L47', 'I48', 'Y49', 'W50', 'A51', 'S52', 'T53', 'R54', 'H55', 'I56', 'G57', 'V58', 'P59', 'D60', 'F62', 'A63', 'G64', 'S65', 'G66', 'S67', 'G68', 'Q89', 'Q90', 'H91', 'Y92', 'S93', 'P94', 'P95', 'W96', 'T97']
pdb_2NZ9 = ['E1', 'I2', 'V3', 'A25', 'S26', 'E27', 'S28', 'V29', 'D30', 'S31', 'Y32', 'G33', 'H34', 'S35', 'F36', 'M37', 'Q38', 'Y53', 'R54', 'A55', 'S56', 'G68', 'S69', 'G70', 'S71', 'G72', 'T73', 'D74', 'F75', 'Q93', 'Q94', 'G95', 'N96', 'E97', 'V98', 'P99', 'F100', 'T101',
            'F27', 'T28', 'F29', 'S30', 'D31', 'H32', 'Y33', 'M34', 'Y35', 'W47', 'A49', 'T50', 'I51', 'S52', 'D53', 'G54', 'G55', 'S56', 'Y57', 'T58', 'Y59', 'Y60', 'S61', 'D62', 'E65', 'T70', 'S71', 'R72', 'D73', 'N74', 'R98', 'Y99', 'R100', 'Y101', 'D102', 'D103', 'A104', 'M105']

pdb_dict = {
    "1BJ1": pdb_1BJ1,
    "1NMB": pdb_1NMB,
    "3BDY": pdb_3BDY,
    "1C08": pdb_1C08,
    "1VFB": pdb_1VFB,
    "3BE1": pdb_3BE1,
    "1CZ8": pdb_1CZ8,
    "1XGP": pdb_1XGP,
    "3BN9": pdb_3BN9,
    "1DQJ": pdb_1DQJ,
    "1XGQ": pdb_1XGQ,
    "3G6D": pdb_3G6D,
    "1DVF": pdb_1DVF,
    "1XGR": pdb_1XGR,
    "3GBN": pdb_3GBN,
    "1JRH": pdb_1JRH,
    "1XGT": pdb_1XGT,
    "3HFM": pdb_3HFM,
    "1KIP": pdb_1KIP,
    "1XGU": pdb_1XGU,
    "3L5X": pdb_3L5X,
    "1KIQ": pdb_1KIQ,
    "1YY9": pdb_1YY9,
    "3N85": pdb_3N85,
    "1KIR": pdb_1KIR,
    "2B2X": pdb_2B2X,
    "3NGB": pdb_3NGB,
    "1MHP": pdb_1MHP,
    "2BDN": pdb_2BDN,
    "3SE8": pdb_3SE8,
    "1MLC": pdb_1MLC,
    "2NY7": pdb_2NY7,
    "4FQY": pdb_4FQY,
    "1N8Z": pdb_1N8Z,
    "2NYY": pdb_2NYY,
    "5C6T": pdb_5C6T,
    "1NCA": pdb_1NCA,
    "2NZ9": pdb_2NZ9
}

data = data[data["LD"] <= 4]
data.reindex()

interfaces = []

for index, row in data.iterrows():
    mutations = re.split(";", row["Mutations"])
    for mutation in mutations:
        mutation = re.sub(r"\w:(\w\d+)\w+", r"\1", mutation)
        is_interface = mutation in pdb_dict[row["#PDB"]]
        if is_interface:
            break
    interfaces.append(is_interface)

data["Interface?"] = interfaces
print(len(data[data["Interface?"] == True]))

try:
    data.to_csv('./raw_datasets/interface_data_use.csv', index=False)
    print("Wrote file.")
except:
    print("Did not output.")