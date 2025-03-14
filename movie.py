mset 1 x180 #由目标分子复制生成180帧（frame）相同的状态（state）;
util.mroll 1, 180, 1 #以当前视野的Y轴旋转分子，得到180个state;
OR util.mrock(1,60,10,1,1)  # issues mdo commands to create +/- 10 deg. rock over 60 frames
set ray_shadows, off  #关掉立体光线阴影
Movie >> 勾选Ray Trace Frames; OR set ray_trace_frames=1
File >> Save Movie As >> PNG Images，选择路径，输入name，会自动保存180个分子旋转过程中的图片，png格式;
OR mpng E:/path/name