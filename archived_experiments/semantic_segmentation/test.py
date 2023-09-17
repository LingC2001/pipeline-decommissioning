import cv2

img = cv2.imread("p14.png",0)
cv2.imshow("",img)
pxl = [0] * 255
for i in range(img.shape[0]):
    for j in range(img.shape[1]):
        pxl[img[i][j]] += 1
print(pxl)

colours = []
for i in range(len(pxl)):
    if pxl[i]!= 0 and i not in colours:
        colours.append(i)
print(colours)

cv2.waitKey(0)