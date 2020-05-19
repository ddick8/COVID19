%% Import data from spreadsheet
% Script for importing data from the following spreadsheet:
%
%    Workbook: D:\University of Toronto\Xiaolin Wei - Covid modeling study\COV_19_Data_revised\Italy_COVID19 updated.xlsx
%    Worksheet: Sheet1
%
% Auto-generated by MATLAB on 14-May-2020 13:33:26

%% Setup the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 10);

% Specify sheet and range
opts.Sheet = "Sheet1";
opts.DataRange = "A4:J81";

% Specify column names and types
opts.VariableNames = ["Var1", "Var2", "Var3", "Var4", "Var5", "cum_confirm_symptomatic", "Var7", "Var8", "Var9", "c_testing"];
opts.SelectedVariableNames = ["cum_confirm_symptomatic", "c_testing"];
opts.VariableTypes = ["char", "char", "char", "char", "char", "double", "char", "char", "char", "double"];

% Specify variable properties
opts = setvaropts(opts, ["Var1", "Var2", "Var3", "Var4", "Var5", "Var7", "Var8", "Var9"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Var1", "Var2", "Var3", "Var4", "Var5", "Var7", "Var8", "Var9"], "EmptyFieldRule", "auto");

% Import the data
tbl = readtable("D:\University of Toronto\Xiaolin Wei - Covid modeling study\COV_19_Data_revised\Italy_COVID19 updated.xlsx", opts, "UseExcel", false);

%% Convert to output type
exp_C = tbl.cum_confirm_symptomatic.';
c_testing = tbl.c_testing.';

%% Clear temporary variables
clear opts tbl