#include <stdio.h>

void print( int arr[], int len )
{
	int i;
	for( i=0; i<len; i++ )
		printf("%i ",arr[i]);
	printf("\n");
}

void Swap( int arr[], int idx1, int idx2 )
{
	int temp = arr[idx1];
	arr[idx1] = arr[idx2];
	arr[idx2] = temp;
}

int Partition( int arr[], int left, int right )
{
	int pivot = arr[left];
	int low = left+1;
	int high = right;

	while( low<=high )
	{
		while( pivot>arr[low] )
			low++;
		while( pivot<arr[high] )
			high--;
		if( low<=high )
			Swap( arr, low, high );
	}

	Swap( arr, left, high );
	return high;
}

void QuickSort( int arr[], int left, int right )
{
	int pivot;

	if( left<=right )
	{
		pivot = Partition( arr, left, right );
		QuickSort( arr, left, pivot-1 );
		QuickSort( arr, pivot+1, right );
	}
}

int main()
{
	int arr[] = { 3, 4, 1, 2, 7 };
	
	print( arr, 5 );

	QuickSort( arr, 0, 4 );

	print( arr, 5 );

	return 0;
}