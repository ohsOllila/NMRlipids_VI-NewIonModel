for f in *.eps
do 
  convert -density 200 -alpha background $f ${f%.eps}.jpg
  epstopdf $f ${f%.eps}-eps-converted-to.pdf
done
