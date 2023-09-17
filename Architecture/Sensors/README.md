# Sensors and other field information:


## Simulators
> We have generated a set of sensor simulators to generate data most similar to those of our in-field sensors. 

### 1. IR scanner

Outputs scans in form of a large matrix. The matrix then is compressed and prepared to be sent through the communication channel.

> 1.1 Compression:
1. gzip : Slow algorithm, high compression
2. numpy "savez_compressed()" : Fast, high compression 


