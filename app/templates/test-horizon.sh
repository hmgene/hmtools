
echo "sample,quality,0,5,10,20,30
A,x,1,2,3,4,5
B,x,5,4,3,2,1" > test.csv 

echo "$data" |   
cat horizon.html | sed s/{{file}}/\"test.csv\"/  > tmp.html
open -a "Google Chrome.app" tmp.html
