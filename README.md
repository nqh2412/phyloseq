# phyloseq
#Hướng dẫn sử dụng DADA2
Link của nhà phát triển package dada2 https://benjjneb.github.io/dada2/tutorial.html


Để xử lý với những dữ liệu sinh học, R-Project là một sự lựa chọn khá tất yếu khi R mang lại một môi trường tính toán thống kê chính xác và linh hoạt. Trong đó Bioconductor là một kho các package của R chuyên để sử lí các số liệu sinh học đặc biệt là dữ liệu trình tự nucleotide. 
Trong nghiên cứu này package Dada2 và phyloseq được áp dụng:

###Mục đích của pipeline
Trong bài hướng dẫn này, OSD-2014 dataset được sử lý qua pipeline DADA2. Bắt đầu từ các trình tự đã được đọc bằng công nghệ đọc trình tự của Ilumina lưu dưới dạng file fastq. Sau khi sử lý qua bằng DADA2, các trình tự của từng mẫu trên được tổng hợp vào bảng ASV - amplicon sequence variant, bảng ASV chứa thông tin về số lượng OTU, và số lần OTU đó có mặt trong mẫu phân tích.

Đồng thời trong bài hướng dẫn này cũng trình bày cách đưa các trình tự này vào pipeline phyloseq và phân tícn microbiome data.

Sau khi cài các package bằng function bioClite(), áp dụng chúng vào R session bằng function library()

###Cài đặt các packages
```{r}
source("https://bioconductor.org/biocLite.R")
biocLite()
biocLite("dada2")
biocLite('phyloseq')
biocLite("DECIPHER")
biocLite("phangorn")
biocLite("BiocStyle")
biocLite("sequencing")
biocLite("BiocGenerics")
biocLite("msa")
```
