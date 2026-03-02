# FastFilter2 Known Bugs

This file lists known issues, their causes, and any recommended workarounds for FastFilter2.

---

## **1. Progress Bar Display Issue with `tqdm`**

**Description:**  
When running FastFilter2 on certain terminals or remote SSH sessions, the stacked `tqdm` progress bars may not render correctly. Symptoms include:

- Overlapping or flickering progress bars  
- Missing updates on some threads  
- Terminal output becomes unreadable during multi-threaded processing  

**Cause:**  
`Tqdm` relies on terminal control sequences to update progress bars. Some terminal environments (e.g., SSH sessions without proper terminal handling) cannot support multi-line or multi-threaded `tqdm` updates.

**Workaround / Fix:**  
- Run FastFilter2 within a `screen` or `tmux` session. This ensures that terminal control sequences are handled correctly and that progress bars render properly.  
  Example:

  ```bash
  screen -S fastfilter
  python fastfilter_pe.py -i /data/project/cutadapt -o /data/project/fastfilter -l 30 -s 30 -p 25 -j 4
