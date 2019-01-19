To execute

1) place the data at "Data"
	a) Binary Data to placed as Data/Binary
	b) Non-Binary Data to placed as Data/NonBinary
2) cd Code
3) ./execute.sh
4) Use the user interface to execute the experiments
	a) In case of FPS you have to provide the number of bits in the fingerprint as input
	   eg : 2048
		this will be considered as target file "targets_2048.fps" and query file "queries_2048.fps"
	b) In the case of Binary and Non-Binary, you will be required to give the file name
	   eg : "DUD"
		this will be considered as target file "DUD" and queries file "DUD_test"
5) A sample Data directory can be downloaded from "https://tinyurl.com/yb5cs2kj"
The results are placed at "Data/Results"


Data format

1) FPS format : The file format used by chemfp, kindly check  http://chemfp.com/fps_format/

2) Binary
	a) each line in the file corresponds a different fingerprint
	b) Data format "{feature_id:1 , feature_id:1 , feature_id:1 , ......}"
	c) E.g.: {1:1 , 2:1 , 3:1 , 4:1 , 5:1 , 6:1 , 7:1 , 8:1 , 9:1 , 10:1 , 11:1 , 12:1 , 13:1 , 14:1 , 15:1 , 16:1 , 17:1 , 18:1 , 19:1 , 20:1 , 21:1 , 22:1 , 23:1 , 24:1 , 25:1}
	d) A fingerprint is identified as the line number in the target file

2) NonBinary
	a) each line in the file corresponds a different fingerprint
	b) Data format "#molecule_id feature_id:feature_value feature_id:feature_value ....."
	c) E.g.: #ZINC00002216 0:4 1:5 2:1 3:1 4:1 5:3 6:1 7:1 8:54 9:23 10:24 11:10 12:12 13:9 14:9 15:16 16:5 17:2 18:4 19:2 20:7 21:2 22:1 23:1 24:2 25:3 26:1 27:7 28:2 29:2 30:3 31:9 32:9 33:4 34:7 35:3 36:1 37:4 38:18 39:3 40:2 41:1 42:1 43:1 44:1 45:5 46:1 47:1 48:77 49:43 50:63 51:20 52:26 53:24 54:4 55:15 56:16 57:34 58:8 59:7 60:19 61:17 62:6 63:6 64:4 65:7 66:18 67:15 68:10 69:4 70:22


