# need  _sbatch_template.sh sbatch_submit.py
# split -l 8 config -d -a 3
for i in x000 x001 x002 x003 x004 x005 x006 x007 x008 x009 x010 x011 x012 x013 x014 x015 x016 x017 x018 x019 x020 x021 x022 x023 x024 x025 x026 x027 x028 x029 x030 x031 x032 x033 x034 x035 x036 x037 x038 x039 x040 x041 x042 x043 x044 x045 x046 x047 x048 x049 x050 x051 x052 x053 x054 x055 x056 x057 x058 x059 x060 x061 x062 x063 x064 x065 x066 x067 x068 x069 x070 x071 x072 x073 x074 x075 x076 x077 x078 x079 x080 x081 x082 x083 x084 x085 x086 x087 x088 x089 x090 x091 x092 x093 x094 x095 x096 x097 x098 x099 x100 x101 x102 x103 x104 x105 x106 x107 x108 x109 x110 x111 x112 x113 x114 x115 x116 x117 x118
do 
dftb_command="/bin/parallel -j 2 -a $i "
echo $dftb_command
python sbatch_submit.py -i "$dftb_command" -j "$i" -f "$i"
sbatch  run_$i.sh
sleep 1
done

