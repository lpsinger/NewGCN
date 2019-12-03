
from django.shortcuts import render

# Create your tests here.

def newgcn(request):
	return render(request, 'table/index.html')
