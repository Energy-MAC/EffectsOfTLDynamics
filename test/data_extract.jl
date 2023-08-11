using Tesseract, Images
image_path = "ss.png"

# Load an image
cd(@__DIR__)
img = Images.load(image_path)

# Perform OCR using Tesseract
extracted_text = ocr(img)

println(extracted_text)
