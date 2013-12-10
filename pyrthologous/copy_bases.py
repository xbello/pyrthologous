import os
import shutil

src_path = "/opt/space/data/Parasitos/secuencias/scripts/to_send/"
tgt_path = src_path + "only_bases/"

for comp in os.listdir(src_path):
    if os.path.isdir(os.path.join(src_path, comp)):
        for fas in os.listdir(os.path.join(src_path, comp)):
            if "bases" in fas:
                if not os.path.isdir(os.path.join(tgt_path, comp)):
                    os.mkdir(os.path.join(tgt_path, comp))
                shutil.copy(os.path.join(src_path, comp, fas),
                    os.path.join(tgt_path, comp, fas))
