\ Model matlab
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 800 C0 - 700 C1 - 600 C2 - 500 C3 - 400 C4 + 500 C5 + 600 C6 + 700 C7
   + 800 C8 + 900 C9
Subject To
 R0: C0 + C1 + C2 + C3 + C4 - C5 - C6 - C7 - C8 - C9 = 0
 R1: - C5 <= -200
 R2: C5 <= 200
 R3: - C6 <= -400
 R4: C6 <= 200
 R5: - C7 <= -600
 R6: C7 <= 200
 R7: - C8 <= -800
 R8: C8 <= 200
 R9: - C9 <= -1000
 R10: C9 <= 200
 R11: - C5 - C6 - C7 - C8 - C9 <= 0
 R12: C5 + C6 + C7 + C8 + C9 <= 4200
 R13: - C0 <= -100
 R14: C0 <= 100
 R15: - C1 <= -200
 R16: C1 <= 100
 R17: - C2 <= -300
 R18: C2 <= 100
 R19: - C3 <= -400
 R20: C3 <= 100
 R21: - C4 <= -500
 R22: C4 <= 100
 R23: - C0 - C1 - C2 - C3 - C4 <= 0
 R24: C0 + C1 + C2 + C3 + C4 <= 2100
Bounds
 C0 free
 C1 free
 C2 free
 C3 free
 C4 free
 C5 free
 C6 free
 C7 free
 C8 free
 C9 free
End
