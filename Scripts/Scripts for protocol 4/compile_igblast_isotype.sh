hook=$1

ls | grep $hook"_" | grep igblast | grep isotype | grep -v all | xargs cat | grep -v Header  > $hook"_igblast_isotype_all.txt"
