tag=$1; e=$2; suf=$3; s=$4; drop=$5; tm=$6; cr=$7; sg=$8
[ "$suf" = "-" ] && suf=''; [ "$drop" = "-" ] && drop=''
cd /workspace; Q22_TAG=$tag Q22_NP=12000 Q22_DROP=$drop Q22_CROSS=$cr Q22_SG=$sg python3 scripts/q22b_shared.py $e "$suf" $s 2>&1 | tail -1
