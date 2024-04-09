sed -n '/Spline/q;p' original/Li-Si.skf >refit/Li-Si.skf
sed -n '/Spline/q;p' original/Si-Si.skf >refit/Si-Si.skf
sed -n '/Spline/q;p' original/Si-Li.skf >refit/Si-Li.skf
sed -n '/Spline/q;p' original/Li-Li.skf >refit/Li-Li.skf
cat Si-Si.spl >> refit/Si-Si.skf
cat Si-Li.spl >> refit/Si-Li.skf
cat Li-Li.spl >> refit/Li-Li.skf
cat Li-Si.spl >> refit/Li-Si.skf
