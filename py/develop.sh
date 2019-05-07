while inotifywait -qqre modify "./py"; do
  clear
  python -m unittest discover
done
