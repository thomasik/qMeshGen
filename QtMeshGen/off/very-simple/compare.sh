for m in block eight genus3 joint knot1
do
  echo $m
  ../../../../MeshCompare/metro407/metro $m/$m.off $m/$m.off.remeshed.off | grep -v Sampling > $m/cmp-metro.txt
done
