Re_s_FB1 = [];
G_FB1 = [];
for i=1:length(FB1.exp)
  Re_s_FB1 = [Re_s_FB1, FB1.exp(i).Re_s];
  G_FB1 = [G_FB1, FB1.exp(i).G];
end
Re_s_FB1
G_FB1

Re_s_FB2 = zeros(length(FB2.exp(1).Re_s), length(FB2.exp));
G_FB2 = zeros(length(FB2.exp(1).Re_s), length(FB2.exp));
for i=1:length(FB2.exp)
  Re_s_FB2(1:length(FB2.exp(i).Re_s), i) = FB2.exp(i).Re_s;
  G_FB2(1:length(FB2.exp(i).G), i) = FB2.exp(i).G;
end
Re_s_FB2
G_FB2
